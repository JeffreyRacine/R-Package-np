/* SCRATCH request-only ANN CDF radius terms. Include after the reviewed stable
 * concentration score and np_plot_nn_response_jet.h; no old owner is changed.
 * Coordinate CDFs retain the canonical point owner. Error is absolute radius
 * term quadrature error, propagated by the R factor before admission. */
/* Own explicit scaled integrand prevents accidental coupling to changes in
 * the NN response-jet quadrature's internal scaling convention. */
typedef struct { double alpha,beta,t,tc,tau,mass; int lower,failed; } NPPlotANNCDFScore;
static void np_plot_ann_cdf_scaled_score(double *x,int n,void *data)
{
  NPPlotANNCDFScore *c=(NPPlotANNCDFScore *)data;
  for(int i=0;i<n;++i) {
    const double prob=c->mass*x[i];
    const double z=qbeta(prob,c->alpha,c->beta,c->lower,0);
    const double zc=qbeta(prob,c->beta,c->alpha,!c->lower,0);
    const double score=np_plot_beta_concentration_score_units(c->t,c->tc,z,zc,c->tau);
    const long double scaled=2.0L*(long double)c->tau*(long double)score;
    x[i]=(double)scaled;
    if(!R_FINITE(x[i])){c->failed=1;x[i]=0.0;}
  }
}

static double np_plot_ann_beta_cdf_radius(double q,double y,double h,
  double lb,double ub,int order,double *absolute_error,double *work,int *iwork)
{
  *absolute_error=0.0;
  if(q<=lb||q>=ub)return 0.0;
  double positive=R_NegInf,negative=R_NegInf;
  for(int j=0;j<order/2;++j) {
    np_beta_shape sh;np_beta_pdf_component coefficient;
    if(np_beta_shape_init(y,y,h,lb,ub,j+1,&sh)!=NP_BETA_OK)
      error("np_plot_ann_unavailable: ANN CDF beta component shape is unavailable");
    if(np_beta_pdf_component_from_shape(&sh,&coefficient)!=NP_BETA_OK||
       np_beta_pdf_component_prepare_coefficient(&coefficient,order,j)!=NP_BETA_OK)
      error("invalid canonical ANN CDF beta component coefficient");
    if(sh.concentration==0.0)continue;
    const double a=1.0+sh.target_unit*sh.concentration;
    const double b=1.0+sh.target_complement_unit*sh.concentration;
    const double L=sh.support_length;
    const double lower=pbeta((q-lb)/L,a,b,1,0);
    const double upper=pbeta((ub-q)/L,b,a,1,0);
    NPPlotANNCDFScore context={a,b,sh.target_unit,sh.target_complement_unit,
      sh.concentration,lower<=upper?lower:upper,lower<=upper,0};
    if(!R_FINITE(context.mass)||context.mass<=0.0||context.mass>1.0)
      error("np_plot_ann_unavailable: ANN CDF beta score tail is not representable");
    double lo=0.0,hi=1.0,epsabs=2e-12,epsrel=2e-10,result=0.0,err=0.0;
    int neval=0,ier=0,limit=256,lenw=4*256,last=0;
    Rdqags(np_plot_ann_cdf_scaled_score,&context,&lo,&hi,&epsabs,&epsrel,
      &result,&err,&neval,&ier,&limit,&lenw,&last,iwork,work);
    if(context.failed||ier!=0||!R_FINITE(result)||!R_FINITE(err)||err<0.0)
      error("np_plot_ann_unavailable: ANN CDF beta concentration quadrature failed");
    /* -h*C_h = 2*tau times the lower-tail score integral. Upper-tail
       integration reverses its sign since the total score integral is zero. */
    const double logscale=log(context.mass)+
      coefficient.log_abs_coefficient;
    if(err>0.0) {
      const double e=exp(logscale+log(err));
      if(!R_FINITE(e)||e==0.0)error("np_plot_ann_unavailable: ANN CDF beta quadrature error is not representable");
      *absolute_error+=e;
    }
    if(result!=0.0) {
      const double term=logscale+log(fabs(result));
      int sign=coefficient.coefficient_sign*(context.lower?1:-1)*(result>0.0?1:-1);
      if(!R_FINITE(term))error("np_plot_ann_unavailable: ANN CDF beta radius is not representable");
      if(sign>0)positive=np_beta_log_add_pair(positive,term);
      else negative=np_beta_log_add_pair(negative,term);
    }
  }
  double absolute=R_NegInf;int sign=0;
  if(!R_FINITE(*absolute_error)||
     np_beta_signed_log_absolute(positive,negative,&absolute,&sign)!=NP_BETA_OK)
    error("np_plot_ann_unavailable: ANN CDF beta signed radius accumulation failed");
  return np_plot_nn_restore_beta(absolute,sign);
}

SEXP C_np_plot_ann_cdf_terms(SEXP train,SEXP query,SEXP bandwidth,
  SEXP kernel,SEXP bounds,SEXP max_bytes)
{
  /* Keep returned attributes rooted across subsequent allocating API calls. */
  SEXP td=PROTECT(getAttrib(train,R_DimSymbol));
  SEXP qd=PROTECT(getAttrib(query,R_DimSymbol));
  SEXP hd=PROTECT(getAttrib(bandwidth,R_DimSymbol));
  SEXP bd=PROTECT(getAttrib(bounds,R_DimSymbol));
  if(TYPEOF(train)!=REALSXP||TYPEOF(query)!=REALSXP||TYPEOF(bandwidth)!=REALSXP||
     TYPEOF(bounds)!=REALSXP||TYPEOF(kernel)!=INTSXP||XLENGTH(kernel)!=1||
     TYPEOF(max_bytes)!=REALSXP||XLENGTH(max_bytes)!=1||
     TYPEOF(td)!=INTSXP||TYPEOF(qd)!=INTSXP||TYPEOF(hd)!=INTSXP||TYPEOF(bd)!=INTSXP||
     XLENGTH(td)!=2||XLENGTH(qd)!=2||XLENGTH(hd)!=2||XLENGTH(bd)!=2)
    error("invalid ANN CDF terms arguments");
  const int n=INTEGER(td)[0],p=INTEGER(td)[1],m=INTEGER(qd)[0],ker=INTEGER(kernel)[0];
  const double count=(double)n*m*p,meta=(double)n*p,limit=REAL(max_bytes)[0];
  const double forecast=24.0*count+128.0*(n+(double)m+p)+16384.0;
  if(n<1||p<1||m<1||INTEGER(qd)[1]!=p||INTEGER(hd)[0]!=n||INTEGER(hd)[1]!=p||
     INTEGER(bd)[0]!=2||INTEGER(bd)[1]!=p||
     (ker<0||(ker>8&&(ker<100||ker>103)))||!R_FINITE(limit)||limit<=0||
     XLENGTH(train)!=(R_xlen_t)meta||XLENGTH(bandwidth)!=(R_xlen_t)meta||
     XLENGTH(query)!=(R_xlen_t)((double)m*p)||XLENGTH(bounds)!=(R_xlen_t)(2.0*p))
    error("invalid ANN CDF terms geometry");
  if(!R_FINITE(forecast)||forecast>limit||count>INT_MAX||meta>INT_MAX)
    error("np_plot_ann_unavailable: native ANN CDF terms exceeds memory budget");
  for(int j=0;j<p;++j) {
    const double lb=REAL(bounds)[2*j],ub=REAL(bounds)[2*j+1];
    if(ISNAN(lb)||ISNAN(ub)||lb>=ub||(ker>=100&&(!R_FINITE(lb)||!R_FINITE(ub))))
      error("invalid ANN CDF bounds");
    for(int i=0;i<n;++i) {
      double y=REAL(train)[i+(size_t)n*j],h=REAL(bandwidth)[i+(size_t)n*j];
      if(!R_FINITE(y)||y<lb||y>ub||!R_FINITE(h)||h<=0)
        error("invalid ANN CDF donor or physical radius");
    }
    for(int q=0;q<m;++q)if(!R_FINITE(REAL(query)[q+(size_t)m*j]))
      error("nonfinite ANN CDF query");
  }
  SEXP out=PROTECT(allocVector(VECSXP,4)),names=PROTECT(allocVector(STRSXP,4));
  const char *keys[]={"cdf","negative.h.derivative","radius.absolute.error","forecast.bytes"};
  SEXP dims=PROTECT(allocVector(INTSXP,3));INTEGER(dims)[0]=n;INTEGER(dims)[1]=m;INTEGER(dims)[2]=p;
  for(int z=0;z<3;++z) {SEXP v=PROTECT(allocVector(REALSXP,(R_xlen_t)count));setAttrib(v,R_DimSymbol,dims);SET_VECTOR_ELT(out,z,v);UNPROTECT(1);}
  SET_VECTOR_ELT(out,3,ScalarReal(forecast));
  for(int z=0;z<4;++z)SET_STRING_ELT(names,z,mkChar(keys[z]));
  setAttrib(out,R_NamesSymbol,names);
  double *work=(double *)R_alloc(1024,sizeof(double));int *iwork=(int *)R_alloc(256,sizeof(int));
  for(int j=0;j<p;++j)for(int i=0;i<n;++i) {
    size_t ij=i+(size_t)n*j;
    double y=REAL(train)[ij],h=REAL(bandwidth)[ij],lb=REAL(bounds)[2*j],ub=REAL(bounds)[2*j+1];
    double a=R_FINITE(lb)?(lb-y)/h:R_NegInf,b=R_FINITE(ub)?(ub-y)/h:R_PosInf;
    double invM=1.0,ak=0.0,bk=0.0;int normalizer_seam=0;
    if(ker<100) {
      NPCKerBoundedNorm norm=np_cker_bounded_norm(ker,y,h,lb,ub);
      invM=norm.inverse_mass;
      if(!R_FINITE(invM)||invM<=0.0)error("np_plot_ann_unavailable: ANN CDF donor normalizer unavailable");
      normalizer_seam=ker==8&&((R_FINITE(a)&&fabs(a)==1.0)||(R_FINITE(b)&&fabs(b)==1.0));
      ak=np_plot_nn_zk(ker,a);bk=np_plot_nn_zk(ker,b);
    }
    for(int q=0;q<m;++q) {
      R_CheckUserInterrupt();double z=REAL(query)[q+(size_t)m*j],C,R=0.0,E=0.0;
      if(z<=lb)C=0.0;
      else if(z>=ub)C=1.0;
      else if(ker==8&&z<y-h)C=0.0;
      else if(ker==8&&z>y+h)C=1.0;
      else if(ker>=100) {
        np_beta_status status=NP_BETA_OK;const int order=2*(ker-99);
        C=np_beta_cdf_order(z,y,h,lb,ub,order,&status);
        if(status!=NP_BETA_OK)error("np_plot_ann_unavailable: ANN canonical beta CDF unavailable");
        R=np_plot_ann_beta_cdf_radius(z,y,h,lb,ub,order,&E,work,iwork);
      } else {
        if(normalizer_seam)error("np_plot_ann_unavailable: ANN CDF uniform normalizer is on a support seam");
        C=np_cker_bounded_cdf(ker+30,z,y,h,lb,ub);
        /* The canonical open K at an isolated empirical uniform radius knot
           defines this population/pilot plug-in. There is no density delta
           face term for an integrated CDF. It is not a finite-map q jet. */
        double u=(z-y)/h,uk=ker==8?u*np_nn_rect(z-y,h):np_plot_nn_zk(ker,u);
        R=(uk-ak+C*(ak-bk))*invM;
      }
      if(!R_FINITE(C)||!R_FINITE(R)||!R_FINITE(E)||E<0.0)
        error("np_plot_ann_unavailable: ANN CDF coordinate or radius term is not representable");
      size_t at=i+(size_t)n*(q+(size_t)m*j);
      REAL(VECTOR_ELT(out,0))[at]=C;REAL(VECTOR_ELT(out,1))[at]=R;REAL(VECTOR_ELT(out,2))[at]=E;
    }
  }
  UNPROTECT(7);return out;
}
