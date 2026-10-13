/* Request-only propagation for the new LSQ uncertainty factor. The protected
 * residual owner remains unchanged. All bounds concern exact real arithmetic
 * on the accepted floating inputs and include the returned result's rounding.
 * Included after np_plot_lp_precision_certificate.h in jksum.c. */
static NPPlotInterval np_plot_lsq_iv_div_positive(NPPlotInterval a,
                                                NPPlotInterval b)
{
  if(!(b.lo>0.0)||!R_FINITE(b.hi))
    error("np_plot_precision_unavailable: residual norm cannot be certified positive");
  NPPlotInterval inverse={nextafter(1.0/b.hi,-INFINITY),
                          nextafter(1.0/b.lo,INFINITY)};
  return np_plot_iv_mul(a,inverse);
}
static NPPlotInterval np_plot_lsq_iv_square(NPPlotInterval a)
{
  double lo=a.lo<=0.0&&a.hi>=0.0?0.0:fmin(fabs(a.lo),fabs(a.hi));
  double hi=fmax(fabs(a.lo),fabs(a.hi));
  NPPlotInterval out={lo==0.0?0.0:fmax(0.0,nextafter(lo*lo,-INFINITY)),
                       hi==0.0?0.0:nextafter(hi*hi,INFINITY)};
  return out;
}
static void np_plot_lsq_precision_matrix(SEXP x,int *m,int *n)
{
  SEXP dim=getAttrib(x,R_DimSymbol);
  if(TYPEOF(x)!=REALSXP||TYPEOF(dim)!=INTSXP||XLENGTH(dim)!=2||
     INTEGER(dim)[0]<1||INTEGER(dim)[1]<1||
     XLENGTH(x)!=(R_xlen_t)INTEGER(dim)[0]*INTEGER(dim)[1])
    error("LSQ precision propagation requires double matrices");
  *m=INTEGER(dim)[0];*n=INTEGER(dim)[1];
  for(R_xlen_t i=0;i<XLENGTH(x);++i)
    if(!R_FINITE(REAL(x)[i]))error("LSQ precision propagation requires finite matrices");
}
static void np_plot_lsq_precision_same(SEXP x,int m,int n,int nonnegative)
{
  int mm,nn;np_plot_lsq_precision_matrix(x,&mm,&nn);
  if(m!=mm||n!=nn)error("LSQ precision propagation dimensions differ");
  if(nonnegative)for(R_xlen_t i=0;i<XLENGTH(x);++i)
    if(REAL(x)[i]<0.0)error("LSQ precision propagation error bound is negative");
}
static NPPlotInterval np_plot_lsq_map_interval(double x,double e)
{
  if(e==0.0)return np_plot_iv(x);
  NPPlotInterval out={nextafter(x-e,-INFINITY),nextafter(x+e,INFINITY)};
  return out;
}
static double np_plot_lsq_precision_ratio(int m,int n,int row,
                                         const double *map,const double *err)
{
  double ms=0.0,es=0.0;
  for(int i=0;i<n;++i){const size_t at=(size_t)i*m+row;
    if(!R_FINITE(map[at])||!R_FINITE(err[at])||err[at]<0.0)return R_PosInf;
    if(fabs(map[at])>ms)ms=fabs(map[at]);if(err[at]>es)es=err[at];}
  if(es==0.0)return 0.0;
  if(ms==0.0||!R_FINITE(es))return R_PosInf;
  double sm=0.0,se=0.0;
  for(int i=0;i<n;++i){const size_t at=(size_t)i*m+row;
    double a=fmax(0.0,nextafter(fabs(map[at])/ms,0.0));
    double b=err[at]==0.0?0.0:nextafter(err[at]/es,INFINITY);
    sm=fmax(0.0,nextafter(sm+fmax(0.0,nextafter(a*a,0.0)),0.0));
    se=np_plot_up_add(se,np_plot_up_mul(b,b));}
  if(sm<=0.0)return R_PosInf;
  return np_plot_up_mul(nextafter(es/ms,INFINITY),
    nextafter(sqrt(nextafter(se/sm,INFINITY)),INFINITY));
}

/* residual: (S,ES,q,index) -> interval for exact (I-S)q / ||I-S||.
 * difference: (upper,Eupper,lower,Elower) -> returned H and enclosed error.
 * factor: (M,EM,residual.interval,actual.A) -> elementwise error and row gate.
 * No fit, callback, global state, MPI collective or numerical owner is used. */
SEXP C_np_plot_lsq_precision_propagate(SEXP mode,SEXP args)
{
  if(TYPEOF(mode)!=STRSXP||XLENGTH(mode)!=1||STRING_ELT(mode,0)==NA_STRING||
     TYPEOF(args)!=VECSXP||XLENGTH(args)!=4)
    error("invalid LSQ precision propagation request");
  const char *kind=CHAR(STRING_ELT(mode,0));
  SEXP map=VECTOR_ELT(args,0),err=VECTOR_ELT(args,1);
  int m,n;np_plot_lsq_precision_matrix(map,&m,&n);
  np_plot_lsq_precision_same(err,m,n,1);
  const double forecast=64.0*(double)m*n+32.0*((double)m+n);
  if(!R_FINITE(forecast)||forecast>128.0*1024.0*1024.0)
    error("np_plot_precision_unavailable: propagation exceeds the resident memory budget");
  if(!strcmp(kind,"residual")){
    SEXP q=VECTOR_ELT(args,2),index=VECTOR_ELT(args,3);
    if(TYPEOF(q)!=REALSXP||XLENGTH(q)!=n||TYPEOF(index)!=INTSXP||XLENGTH(index)!=m)
      error("invalid LSQ residual precision vectors");
    for(int i=0;i<n;++i)if(!R_FINITE(REAL(q)[i]))error("nonfinite LSQ precision response");
    for(int j=0;j<m;++j)if(INTEGER(index)[j]<1||INTEGER(index)[j]>n)
      error("invalid LSQ residual precision row index");
    SEXP out=PROTECT(allocMatrix(REALSXP,m,2));
    for(int j=0;j<m;++j){
      double scale=0.0;
      for(int i=0;i<n;++i){const size_t at=(size_t)i*m+j;
        NPPlotInterval e=np_plot_iv_add(np_plot_iv(i==INTEGER(index)[j]-1?1.0:0.0),
          np_plot_iv_neg(np_plot_lsq_map_interval(REAL(map)[at],REAL(err)[at])));
        const double v=np_plot_iv_abs(e);if(v>scale)scale=v;}
      if(!(scale>0.0)||!R_FINITE(scale))
        error("np_plot_precision_unavailable: residual scale cannot be certified");
      NPPlotInterval num=np_plot_iv(0.0),ss=np_plot_iv(0.0);
      for(int i=0;i<n;++i){const size_t at=(size_t)i*m+j;
        NPPlotInterval e=np_plot_iv_add(np_plot_iv(i==INTEGER(index)[j]-1?1.0:0.0),
          np_plot_iv_neg(np_plot_lsq_map_interval(REAL(map)[at],REAL(err)[at])));
        e=np_plot_lsq_iv_div_positive(e,np_plot_iv(scale));
        num=np_plot_iv_add(num,np_plot_iv_mul(e,np_plot_iv(REAL(q)[i])));
        ss=np_plot_iv_add(ss,np_plot_lsq_iv_square(e));}
      if(!(ss.lo>0.0)||!R_FINITE(ss.hi))
        error("np_plot_precision_unavailable: residual norm cannot be certified positive");
      NPPlotInterval norm={nextafter(sqrt(ss.lo),-INFINITY),nextafter(sqrt(ss.hi),INFINITY)};
      NPPlotInterval z=np_plot_lsq_iv_div_positive(num,norm);
      if(!R_FINITE(z.lo)||!R_FINITE(z.hi))
        error("np_plot_precision_unavailable: normalized residual interval is not representable");
      REAL(out)[j]=z.lo;REAL(out)[j+m]=z.hi;
    }
    UNPROTECT(1);return out;
  }
  if(strcmp(kind,"difference")&&strcmp(kind,"factor"))
    error("unknown LSQ precision propagation mode");
  SEXP other=VECTOR_ELT(args,2),last=VECTOR_ELT(args,3);
  if(!strcmp(kind,"difference")){
    np_plot_lsq_precision_same(other,m,n,0);np_plot_lsq_precision_same(last,m,n,1);
  }else{
    np_plot_lsq_precision_same(other,n,2,0);np_plot_lsq_precision_same(last,m,n,0);
    for(int i=0;i<n;++i)if(REAL(other)[i]>REAL(other)[i+n])
      error("invalid LSQ normalized residual interval");
  }
  SEXP value=PROTECT(allocMatrix(REALSXP,m,n)),error_bound=PROTECT(allocMatrix(REALSXP,m,n));
  SEXP ratio=PROTECT(allocVector(REALSXP,m)),certified=PROTECT(allocVector(LGLSXP,m));
  for(int i=0;i<n;++i)for(int j=0;j<m;++j){const size_t at=(size_t)i*m+j;
    NPPlotInterval truth=np_plot_lsq_map_interval(REAL(map)[at],REAL(err)[at]);
    double v;
    if(!strcmp(kind,"difference")){
      truth=np_plot_iv_add(truth,np_plot_iv_neg(np_plot_lsq_map_interval(REAL(other)[at],REAL(last)[at])));
      v=REAL(map)[at]-REAL(other)[at];
    }else{
      NPPlotInterval z={REAL(other)[i],REAL(other)[i+n]};
      truth=np_plot_iv_mul(truth,z);v=REAL(last)[at];
    }
    if(!R_FINITE(truth.lo)||!R_FINITE(truth.hi)||truth.lo>truth.hi||!R_FINITE(v))
      error("np_plot_precision_unavailable: propagated interval is not representable");
    const double bound=np_plot_iv_abs(np_plot_iv_add(truth,np_plot_iv(-v)));
    if(!R_FINITE(bound)||bound<0.0)
      error("np_plot_precision_unavailable: propagated error is not representable");
    REAL(value)[at]=v;REAL(error_bound)[at]=bound;
  }
  for(int j=0;j<m;++j){double r=np_plot_lsq_precision_ratio(m,n,j,REAL(value),REAL(error_bound));
    REAL(ratio)[j]=r;LOGICAL(certified)[j]=R_FINITE(r)&&r<=2e-11;}
  SEXP out=PROTECT(allocVector(VECSXP,4)),names=PROTECT(allocVector(STRSXP,4));
  SEXP values[4]={value,error_bound,ratio,certified};
  const char *keys[]={"value","error","relative.error","certified"};
  for(int i=0;i<4;++i){SET_VECTOR_ELT(out,i,values[i]);SET_STRING_ELT(names,i,mkChar(keys[i]));}
  setAttrib(out,R_NamesSymbol,names);UNPROTECT(6);return out;
}
