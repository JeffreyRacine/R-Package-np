/* SCRATCH REQUEST-ONLY PROTOTYPE. Include after np_plot_nn_response_jet.h
 * and the independently reviewed np_plot_beta_concentration_score.h,
 * then np_plot_beta_pdf_agreement.h.
 * No existing fit, uncertainty, kernel, or selector owner calls this entry.
 * h is the canonical donor-owned ANN physical radius, not an NN rank.
 * A zero-row query requests only bounded mass metadata (moving-face pilot).
 */

/* Each higher-order beta kernel is the canonical signed combination of
 * scale components. Differentiate each concentration, then combine in signed
 * log space. Reuse coefficient preparation; no copied coefficient table. */
static double np_plot_ann_beta_radius_order(double q,double y,double h,
  double lb,double ub,int order,double *weight_error,double *radius_error)
{
  *weight_error=0.0;*radius_error=0.0;
  if(!np_beta_order_supported(order))error("unsupported canonical beta order");
  double positive=R_NegInf,negative=R_NegInf;
  for(int j=0;j<order/2;++j) {
    np_beta_shape s;np_beta_pdf_component component;
    np_beta_status status=np_beta_shape_init(q,y,h,lb,ub,j+1,&s);
    if(status!=NP_BETA_OK)error("np_plot_ann_unavailable: ANN beta radius component shape failed");
    if(np_beta_pdf_component_from_shape(&s,&component)!=NP_BETA_OK||
       np_beta_pdf_component_prepare_coefficient(&component,order,j)!=NP_BETA_OK)
      error("ANN beta radius canonical coefficient failed");
    if(s.concentration==0.0)continue;
    double logterm,score;int sign=component.coefficient_sign;
    if(s.observation_endpoint!=0) {
      /* A fixed endpoint donor contributes identically zero except when
         query and donor are the same endpoint. There f=(1+tau)/L. */
      if(!((y==lb&&q==lb)||(y==ub&&q==ub)))continue;
      logterm=log(2.0)+log(s.concentration)-log(s.support_length);
    } else {
      double lf=np_beta_log_pdf_scale(&s,&status);
      if(status!=NP_BETA_OK)error("np_plot_ann_unavailable: ANN beta radius component PDF failed");
      double discrepancy=np_plot_beta_pdf_component_discrepancy(&s,lf);
      double coefficient_abs=exp(component.log_abs_coefficient);
      double weighted_discrepancy=coefficient_abs*discrepancy;
      if(!R_FINITE(weighted_discrepancy)||(discrepancy>0.0&&weighted_discrepancy==0.0))
        error("np_plot_beta_numerical_unavailable: component discrepancy scaling failed");
      *weight_error+=weighted_discrepancy;
      if(lf==R_NegInf) {
        if(discrepancy!=0.0)error("np_plot_beta_numerical_unavailable: canonical component vanished");
        continue;
      }
      score=np_plot_beta_concentration_score(q,y,lb,ub,s.concentration);
      if(!R_FINITE(score))error("np_plot_ann_unavailable: ANN beta radius score is not representable");
      const long double radius_multiplier=2.0L*(long double)s.concentration*fabsl((long double)score);
      const long double radius_discrepancy=radius_multiplier*(long double)weighted_discrepancy;
      const double radius_discrepancy_double=(double)radius_discrepancy;
      if(radius_discrepancy>0.0L&&radius_discrepancy_double==0.0)
        error("np_plot_beta_numerical_unavailable: radius discrepancy underflow");
      *radius_error+=radius_discrepancy_double;
      if(!R_FINITE(*weight_error)||!R_FINITE(*radius_error))
        error("np_plot_beta_numerical_unavailable: signed component discrepancy overflow");
      if(score==0.0)continue;
      if(score<0.0)sign=-sign;
      logterm=lf+log(2.0)+log(s.concentration)+log(fabs(score));
    }
    logterm+=component.log_abs_coefficient;
    if(!R_FINITE(logterm))error("np_plot_ann_unavailable: ANN beta radius component is not representable");
    if(sign>0)positive=np_beta_log_add_pair(positive,logterm);
    else negative=np_beta_log_add_pair(negative,logterm);
  }
  double absolute=R_NegInf;int sign=0;
  if(np_beta_signed_log_absolute(positive,negative,&absolute,&sign)!=NP_BETA_OK)
    error("np_plot_ann_unavailable: ANN beta radius signed accumulation failed");
  return np_plot_nn_restore_beta(absolute,sign);
}

SEXP C_np_plot_ann_kernel_terms(SEXP train,SEXP query,SEXP bandwidth,
  SEXP kernel,SEXP bounds,SEXP max_bytes)
{
  SEXP td=getAttrib(train,R_DimSymbol),qd=getAttrib(query,R_DimSymbol),
    hd=getAttrib(bandwidth,R_DimSymbol),bd=getAttrib(bounds,R_DimSymbol);
  if(TYPEOF(train)!=REALSXP||TYPEOF(query)!=REALSXP||TYPEOF(bandwidth)!=REALSXP||
     TYPEOF(bounds)!=REALSXP||TYPEOF(kernel)!=INTSXP||XLENGTH(kernel)!=1||
     TYPEOF(max_bytes)!=REALSXP||XLENGTH(max_bytes)!=1||
     TYPEOF(td)!=INTSXP||TYPEOF(qd)!=INTSXP||TYPEOF(hd)!=INTSXP||TYPEOF(bd)!=INTSXP||
     XLENGTH(td)!=2||XLENGTH(qd)!=2||XLENGTH(hd)!=2||XLENGTH(bd)!=2)
    error("invalid ANN plot kernel terms arguments");
  const int n=INTEGER(td)[0],p=INTEGER(td)[1],m=INTEGER(qd)[0],ker=INTEGER(kernel)[0];
  const double count=(double)n*m*p,meta=(double)n*p,limit=REAL(max_bytes)[0];
  const double forecast=32.0*count+16.0*meta+128.0*(n+(double)m+p)+4096.0;
  if(n<1||p<1||m<0||INTEGER(qd)[1]!=p||INTEGER(hd)[0]!=n||INTEGER(hd)[1]!=p||
     INTEGER(bd)[0]!=2||INTEGER(bd)[1]!=p||
     (ker<0||(ker>8&&(ker<100||ker>103)))||!R_FINITE(limit)||limit<=0||
     XLENGTH(train)!=(R_xlen_t)meta||XLENGTH(bandwidth)!=(R_xlen_t)meta||
     XLENGTH(query)!=(R_xlen_t)((double)m*p)||XLENGTH(bounds)!=(R_xlen_t)(2.0*p))
    error("invalid ANN plot kernel terms geometry");
  if(!R_FINITE(forecast)||forecast>limit||count>INT_MAX||meta>INT_MAX)
    error("np_plot_ann_unavailable: native ANN plot kernel terms exceeds memory budget");
  for(int j=0;j<p;++j) {
    const double lb=REAL(bounds)[2*j],ub=REAL(bounds)[2*j+1];
    if(ISNAN(lb)||ISNAN(ub)||lb>=ub||(ker>=100&&(!R_FINITE(lb)||!R_FINITE(ub))))
      error("invalid ANN plot kernel bounds");
    for(int i=0;i<n;++i) {
      double y=REAL(train)[i+(size_t)n*j],h=REAL(bandwidth)[i+(size_t)n*j];
      if(!R_FINITE(y)||y<lb||y>ub||!R_FINITE(h)||h<=0)
        error("invalid ANN plot donor or physical radius");
    }
    for(int q=0;q<m;++q) {
      double z=REAL(query)[q+(size_t)m*j];
      if(!R_FINITE(z)||z<lb||z>ub)error("invalid ANN plot query");
      if(ker==8&&(z<=lb||z>=ub))error("np_plot_ann_unavailable: uniform query coincides with physical support face");
    }
  }
  SEXP out=PROTECT(allocVector(VECSXP,7)),names=PROTECT(allocVector(STRSXP,7));
  const char *keys[]={"weight","negative.h.derivative","mass","mass.log.h.derivative","weight.absolute.discrepancy","radius.absolute.discrepancy","forecast.bytes"};
  SEXP dims=PROTECT(allocVector(INTSXP,3));INTEGER(dims)[0]=n;INTEGER(dims)[1]=m;INTEGER(dims)[2]=p;
  for(int z=0;z<2;++z) {SEXP x=PROTECT(allocVector(REALSXP,(R_xlen_t)count));setAttrib(x,R_DimSymbol,dims);SET_VECTOR_ELT(out,z,x);UNPROTECT(1);}
  for(int z=2;z<4;++z)SET_VECTOR_ELT(out,z,allocMatrix(REALSXP,n,p));
  for(int z=4;z<6;++z) {SEXP v=PROTECT(allocVector(REALSXP,(R_xlen_t)count));setAttrib(v,R_DimSymbol,dims);SET_VECTOR_ELT(out,z,v);UNPROTECT(1);}
  SET_VECTOR_ELT(out,6,ScalarReal(forecast));
  for(int z=0;z<7;++z)SET_STRING_ELT(names,z,mkChar(keys[z]));setAttrib(out,R_NamesSymbol,names);
  for(int j=0;j<p;++j)for(int i=0;i<n;++i) {
    R_CheckUserInterrupt();size_t ij=i+(size_t)n*j;
    double y=REAL(train)[ij],h=REAL(bandwidth)[ij],lb=REAL(bounds)[2*j],ub=REAL(bounds)[2*j+1];
    double a=R_FINITE(lb)?(lb-y)/h:R_NegInf,b=R_FINITE(ub)?(ub-y)/h:R_PosInf;
    double invM=1.0,scale=1.0/h,c=0.0;
    if(ker<100) {
      NPCKerBoundedNorm norm=np_cker_bounded_norm(ker,y,h,lb,ub);
      invM=norm.inverse_mass;scale=norm.density_scale;
      if(!R_FINITE(invM)||invM<=0.0||!R_FINITE(scale)||scale<=0.0)
        error("np_plot_ann_unavailable: ANN plot donor normalizer is unavailable");
      /* The uniform normalizer has a kink when its moving support meets a
         fixed physical bound. Smooth kernels use their canonical derivatives. */
      if(ker==8&&((R_FINITE(a)&&fabs(a)==1.0)||(R_FINITE(b)&&fabs(b)==1.0)))
        error("np_plot_ann_unavailable: ANN plot uniform normalizer is on a support seam");
      c=(np_plot_nn_zk(ker,a)-np_plot_nn_zk(ker,b))*invM;
    }
    REAL(VECTOR_ELT(out,2))[ij]=1.0/invM;
    REAL(VECTOR_ELT(out,3))[ij]=c;
    for(int q=0;q<m;++q) {
      double z=REAL(query)[q+(size_t)m*j],w,r,we=0.0,re=0.0;
      if(ker>=100) {
        np_beta_status status=NP_BETA_OK;int sign=0;
        const int beta_order=2*(ker-99);
        double lf=np_beta_log_abs_pdf_order(z,y,h,lb,ub,beta_order,&sign,&status);
        if(status!=NP_BETA_OK)error("np_plot_ann_unavailable: ANN plot beta PDF is unavailable");
        w=np_plot_nn_restore_beta(lf,sign);
        r=np_plot_ann_beta_radius_order(z,y,h,lb,ub,beta_order,&we,&re);
      } else {
        double u=(z-y)/h,ku=0.0;
        /* Isolated empirical query faces retain the canonical open-support
           zero. The aggregate distributional face term is assembled by the
           ANN density influence owner; this is not a fitted-map q derivative. */
        double base=ker==8?np_nn_rect(z-y,h):np_plot_nn_k(ker,u);
        w=base*scale;
        if(ker!=8)np_ckernelv(ker+20,&y,1,0,z,h,&ku,NULL,0,1,1,1.0,0,R_NegInf,R_PosInf,NULL,NULL);
        r=w+(z-y)*ku*scale/h+w*c;
      }
      if(!R_FINITE(w)||!R_FINITE(r))error("np_plot_ann_unavailable: ANN plot kernel radius term is not representable");
      size_t at=i+(size_t)n*(q+(size_t)m*j);
      REAL(VECTOR_ELT(out,0))[at]=w;REAL(VECTOR_ELT(out,1))[at]=r;
      REAL(VECTOR_ELT(out,4))[at]=we;REAL(VECTOR_ELT(out,5))[at]=re;
    }
  }
  UNPROTECT(3);return out;
}
