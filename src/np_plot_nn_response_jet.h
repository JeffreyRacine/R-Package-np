/* Scratch, request-only response geometry jet. Include after jksum.c definitions.
 * No established owner calls this entry. Public qreg response order is two. */
#include <R_ext/Applic.h>

/* New-request absolute discrepancy arithmetic. Preserve nonzero uncertainty
 * rather than losing it before a later derivative amplification. */
static double np_plot_nn_error_product(double a,double b)
{
  if(a==0.0||b==0.0)return 0.0;
  double z=a*b;
  if(!R_FINITE(z)||z<=0.0)
    error("np_plot_beta_numerical_unavailable: response discrepancy multiplication is not representable");
  z=nextafter(z,INFINITY);
  if(!R_FINITE(z))error("np_plot_beta_numerical_unavailable: response discrepancy overflow");
  return z;
}
static double np_plot_nn_error_sum(double a,double b)
{
  if(a==0.0)return b;if(b==0.0)return a;
  double z=nextafter(a+b,INFINITY);
  if(!R_FINITE(z)||z<0.0)error("np_plot_beta_numerical_unavailable: response discrepancy sum is not representable");
  return z;
}
static double np_plot_nn_log_multiplier_error(double e,double log_multiplier)
{
  if(e==0.0||log_multiplier==R_NegInf)return 0.0;
  double v=exp(log(e)+log_multiplier);
  if(!R_FINITE(v)||v==0.0)
    error("np_plot_beta_numerical_unavailable: response derivative discrepancy is not representable");
  v=nextafter(v,INFINITY);
  if(!R_FINITE(v))error("np_plot_beta_numerical_unavailable: response derivative discrepancy overflow");
  return v;
}

typedef struct {
  double alpha, beta, t, tc, tau, mass;
  int lower, failed;
} NPPlotNNBetaScore;

static void np_plot_nn_beta_score_integrand(double *x, int n, void *data)
{
  NPPlotNNBetaScore *c=(NPPlotNNBetaScore *)data;
  for(int i=0;i<n;++i) {
    const double p=c->mass*x[i];
    /* Work with both tails directly; do not subtract a rounded quantile from 1. */
    const double z=qbeta(p,c->alpha,c->beta,c->lower,0);
    const double zc=qbeta(p,c->beta,c->alpha,!c->lower,0);
    x[i]=2.0*(c->tau*np_plot_beta_concentration_score_units(c->t,c->tc,z,zc,c->tau));
    if(!R_FINITE(x[i])) {c->failed=1;x[i]=0.0;}
  }
}

/* Integral of the concentration score, not a finite difference in bandwidth.
 * Error returned on the physical C_h scale. Caller retains it in the receipt. */
static double np_plot_nn_beta_cdf_h(double q, double y, double h,
  double lb, double ub, double *absolute_error, double *work, int *iwork)
{
  np_beta_shape sh;
  if(np_beta_shape_init(y,y,h,lb,ub,1,&sh)!=NP_BETA_OK)
    error("np_plot_nn_unavailable: NN response beta concentration is unavailable");
  if(sh.concentration==0.0) { *absolute_error=0.0; return 0.0; }
  const double a=1.0+sh.target_unit*sh.concentration;
  const double b=1.0+sh.target_complement_unit*sh.concentration;
  const double L=sh.support_length;
  const double lower=pbeta((q-lb)/L,a,b,1,0);
  const double upper=pbeta((ub-q)/L,b,a,1,0);
  NPPlotNNBetaScore context={a,b,sh.target_unit,sh.target_complement_unit,
    sh.concentration,
    lower<=upper ? lower : upper,lower<=upper,0};
  if(!R_FINITE(context.mass)||context.mass<0.0||context.mass>1.0)
    error("np_plot_nn_unavailable: NN response beta tail probability is unavailable");
  if(context.mass==0.0) error("np_plot_nn_unavailable: NN response beta score tail is not representable");
  double lo=0.0,hi=1.0,epsabs=2e-12,epsrel=2e-10,result=0.0,err=0.0;
  int neval=0,ier=0,limit=256,lenw=4*256,last=0;
  Rdqags(np_plot_nn_beta_score_integrand,&context,&lo,&hi,&epsabs,&epsrel,
    &result,&err,&neval,&ier,&limit,&lenw,&last,iwork,work);
  const double coefficient=-(context.mass/h)*
    (context.lower ? 1.0 : -1.0);
  const double value=coefficient*result;
  if(coefficient==0.0||(value==0.0&&result!=0.0))
    error("np_plot_nn_unavailable: NN beta concentration integral scale is not representable");
  *absolute_error=np_plot_nn_error_product(fabs(coefficient),err);
  if(context.failed||ier!=0||!R_FINITE(value)||!R_FINITE(*absolute_error)||
     *absolute_error>2e-9*(1.0+fabs(value)))
    error("np_plot_nn_unavailable: NN response beta concentration integral did not meet its error bound");
  return value;
}

static double np_plot_nn_beta_pdf_h(double q,double y,double h,double lb,double ub,double f)
{
  np_beta_shape sh;
  if(np_beta_shape_init(q,y,h,lb,ub,1,&sh)!=NP_BETA_OK)
    error("np_plot_nn_unavailable: NN response beta density concentration is unavailable");
  /* With an interior query, an endpoint observation has zero density for
     every nearby h, so the composite derivative is exactly zero. */
  if(f==0.0 && (y==lb||y==ub)) return 0.0;
  double score=np_plot_beta_concentration_score(q,y,lb,ub,sh.concentration);
  return f*((-2.0*(sh.concentration*score))/h);
}

static double np_plot_nn_restore_beta(double logabs,int sign)
{
  if(sign==0) return 0.0;
  double out=sign*exp(logabs);
  if(!R_FINITE(out) || (out==0.0 && logabs!=-INFINITY))
    error("np_plot_nn_unavailable: NN response beta value or derivative is not representable");
  return out;
}

static double np_plot_nn_k(int kernel,double z)
{
  if(!R_FINITE(z)) return 0.0;
  return np_cker_base_eval(kernel,z);
}

static double np_plot_nn_zk(int kernel,double z)
{
  return R_FINITE(z) ? z*np_plot_nn_k(kernel,z) : 0.0;
}

/* Compact support faces are regular only when the composed path is tangent.
 * This admits the persistent excluded donor on a uniform GNN radius face. */
static void np_plot_nn_face(int kernel,double z,double velocity)
{
  if(kernel==4||kernel==8) {
    const double edge=cksup[kernel][1];
    if(fabs(z)==edge && velocity!=0.0)
      error("np_plot_nn_unavailable: NN response derivative is on a nonregular compact support face");
  }
}

static double np_plot_nn_radius_slope(const double *y,int n,double q,
  double radius,int k,double scale)
{
  int below=0,left=0,right=0;
  for(int i=0;i<n;++i) {
    double d=fabs(q-y[i]);
    if(d<radius) ++below;
    else if(d==radius) { if(q>y[i])++left;else if(q<y[i])++right; }
  }
  int r=k-below;
  if(r<1||r>left+right)
    error("NN response radius did not match canonical observation geometry");
  int slope_right=r<=right ? -1 : 1;
  int slope_left=r<=left ? 1 : -1;
  if(slope_right!=slope_left)
    error("np_plot_nn_unavailable: NN response radius has different left and right query derivatives");
  return scale*slope_right;
}

SEXP C_np_plot_nn_response_jet(SEXP train,SEXP query,SEXP requested,
  SEXP mode,SEXP kernel,SEXP bounds,SEXP max_bytes)
{
  if(TYPEOF(train)!=REALSXP||TYPEOF(query)!=REALSXP||TYPEOF(requested)!=REALSXP||
     XLENGTH(requested)!=1||TYPEOF(mode)!=INTSXP||XLENGTH(mode)!=1||
     TYPEOF(kernel)!=INTSXP||XLENGTH(kernel)!=1||TYPEOF(bounds)!=REALSXP||
     XLENGTH(bounds)!=2||TYPEOF(max_bytes)!=REALSXP||XLENGTH(max_bytes)!=1||
     XLENGTH(train)<2||XLENGTH(query)<1||XLENGTH(train)>INT_MAX||XLENGTH(query)>INT_MAX)
    error("invalid NN response jet arguments");
  int n=(int)XLENGTH(train),m=(int)XLENGTH(query),bw=INTEGER(mode)[0],ker=INTEGER(kernel)[0];
  double lb=REAL(bounds)[0],ub=REAL(bounds)[1],limit=REAL(max_bytes)[0];
  double count=(double)n*m,forecast=64.0*count+64.0*((double)n+(double)m)+1024.0*8+256.0*4;
  if((bw!=BW_GEN_NN&&bw!=BW_ADAP_NN)||(ker!=0&&ker!=4&&ker!=8&&ker!=100)||
     ISNAN(lb)||ISNAN(ub)||lb>=ub||!R_FINITE(limit)||limit<=0||
     !R_FINITE(forecast)||forecast>limit||count>INT_MAX||
     (ker==100&&(!R_FINITE(lb)||!R_FINITE(ub))))
    error("invalid or oversized NN response jet geometry");
  for(int i=0;i<n;++i)
    if(!R_FINITE(REAL(train)[i])||REAL(train)[i]<lb||REAL(train)[i]>ub)
      error("NN response donors must be finite and in support");
  for(int j=0;j<m;++j)
    if(!R_FINITE(REAL(query)[j])||REAL(query)[j]<=lb||REAL(query)[j]>=ub)
      error("NN response queries must be strictly inside support");
  int k=0;double scale=0.0;
  if(np_nn_lookup_from_scale(n,1,REAL(requested)[0],&k,&scale,NULL))
    error("invalid NN response count");
  SEXP out=PROTECT(allocVector(VECSXP,10)),names=PROTECT(allocVector(STRSXP,10));
  const char *keys[]={"cdf","pdf","cdf.q","pdf.q","cdf.q.error","bandwidth","bandwidth.q","forecast.bytes","pdf.absolute.discrepancy","pdf.q.absolute.discrepancy"};
  for(int z=0;z<5;++z) SET_VECTOR_ELT(out,z,allocMatrix(REALSXP,m,n));
  SET_VECTOR_ELT(out,5,allocVector(REALSXP,bw==BW_GEN_NN?m:n));
  SET_VECTOR_ELT(out,6,allocVector(REALSXP,m));
  SET_VECTOR_ELT(out,7,ScalarReal(forecast));
  SET_VECTOR_ELT(out,8,allocMatrix(REALSXP,m,n));
  SET_VECTOR_ELT(out,9,allocMatrix(REALSXP,m,n));
  for(int z=0;z<10;++z) SET_STRING_ELT(names,z,mkChar(keys[z]));
  setAttrib(out,R_NamesSymbol,names);
  double *ht=REAL(VECTOR_ELT(out,5)),*hp=REAL(VECTOR_ELT(out,6));
  double *raw_radius=(double *)R_alloc(bw==BW_GEN_NN?m:n,sizeof(double));
  NPNNGeometryContext geometry={.mode=NP_NN_QUERY_EXTERNAL,.eval_to_train=NULL,.adaptive_successor=NULL};
  NPNNGeometryStatus status=NP_NN_GEOMETRY_OK;
  status=bw==BW_GEN_NN ?
    compute_nn_distance_train_eval_ctx(n,m,1,REAL(train),REAL(query),k,&geometry,raw_radius) :
    compute_nn_distance_status(n,1,REAL(train),k,raw_radius);
  if(status==NP_NN_GEOMETRY_ZERO_RADIUS || status==NP_NN_GEOMETRY_NONFINITE_RADIUS)
    error("np_plot_nn_unavailable: NN response canonical radius is not representable (%d)",(int)status);
  if(status!=NP_NN_GEOMETRY_OK)
    error("NN response canonical radius preparation failed (%d)",(int)status);
  for(int z=0;z<(bw==BW_GEN_NN?m:n);++z) {
    ht[z]=scale*raw_radius[z];
    if(!R_FINITE(ht[z])||ht[z]<=0.0)error("np_plot_nn_unavailable: NN response physical radius is unavailable");
  }
  double *work=(double *)R_alloc(1024,sizeof(double));
  int *iwork=(int *)R_alloc(256,sizeof(int));
  for(int j=0;j<m;++j) {
    R_CheckUserInterrupt();
    double q=REAL(query)[j];
    hp[j]=bw==BW_GEN_NN ? np_plot_nn_radius_slope(REAL(train),n,q,raw_radius[j],k,scale) : 0.0;
    for(int i=0;i<n;++i) {
      double y=REAL(train)[i],h=ht[bw==BW_GEN_NN?j:i],v=hp[j],C,f,Cq,fq,error_bound=0.0,pdf_error=0.0,pdfq_error=0.0;
      if(ker==100) {
        np_beta_status bs=NP_BETA_OK;int sign=0;
        C=np_beta_cdf_order(q,y,h,lb,ub,2,&bs);
        if(bs!=NP_BETA_OK) error("NN beta CDF primitive failed");
        double log_f=np_beta_log_abs_pdf_order(q,y,h,lb,ub,2,&sign,&bs);
        if(bs!=NP_BETA_OK) error("NN beta PDF primitive failed");
        f=np_plot_nn_restore_beta(log_f,sign);
        np_beta_shape pdf_shape;
        if(np_beta_shape_init(q,y,h,lb,ub,1,&pdf_shape)!=NP_BETA_OK)
          error("NN beta PDF agreement shape failed");
        pdf_error=np_plot_beta_pdf_component_discrepancy(&pdf_shape,log_f);
        double log_Fq=np_beta_log_abs_pdf_order(y,q,h,lb,ub,2,&sign,&bs);
        if(bs!=NP_BETA_OK) error("NN beta CDF query derivative failed");
        Cq=np_plot_nn_restore_beta(log_Fq,sign);
        np_beta_shape cdf_shape;
        if(np_beta_shape_init(y,q,h,lb,ub,1,&cdf_shape)!=NP_BETA_OK)
          error("NN beta CDF slope agreement shape failed");
        error_bound=np_plot_beta_pdf_component_discrepancy(&cdf_shape,log_Fq);
        np_beta_derivative d;
        bs=np_beta_pdf_derivative_order(q,y,h,lb,ub,2,&d);
        if(bs!=NP_BETA_OK||d.jump_sign!=0) error("np_plot_nn_unavailable: NN beta PDF query derivative is unavailable");
        fq=np_plot_nn_restore_beta(d.regular_log_absolute,d.regular_sign);
        if(pdf_error>0.0 && d.regular_sign!=0)
          pdfq_error=np_plot_nn_log_multiplier_error(pdf_error,d.regular_log_absolute-log_f);
        if(v!=0.0) {
          double quadrature_error=0.0;
          Cq+=v*np_plot_nn_beta_cdf_h(q,y,h,lb,ub,&quadrature_error,work,iwork);
          error_bound=np_plot_nn_error_sum(error_bound,np_plot_nn_error_product(fabs(v),quadrature_error));
          fq+=v*np_plot_nn_beta_pdf_h(q,y,h,lb,ub,f);
          const double score=np_plot_beta_concentration_score(q,y,lb,ub,pdf_shape.concentration);
          const double multiplier=fabs((-2.0*(pdf_shape.concentration*score))/h);
          if(!R_FINITE(multiplier))error("np_plot_beta_numerical_unavailable: response radius multiplier is unavailable");
          pdfq_error=np_plot_nn_error_sum(pdfq_error,np_plot_nn_error_product(fabs(v),
            np_plot_nn_error_product(pdf_error,multiplier)));
        }
      } else {
        double u=(q-y)/h,a=R_FINITE(lb)?(lb-y)/h:R_NegInf,b=R_FINITE(ub)?(ub-y)/h:R_PosInf;
        double M=np_cker_interval_mass(ker,a,b);
        if(!R_FINITE(M)||M<=0.0)error("np_plot_nn_unavailable: NN response donor normalizer is unavailable");
        C=np_cker_bounded_cdf(ker+30,q,y,h,lb,ub);
        np_plot_nn_face(ker,u,1.0-u*v);
        np_plot_nn_face(ker,a,-a*v);np_plot_nn_face(ker,b,-b*v);
        double K=(ker==8)?np_nn_rect(q-y,h):np_plot_nn_k(ker,u);
        double Ku;
        np_ckernelv(ker+20,&y,1,0,q,h,&Ku,NULL,0,1,1,1.0,0,R_NegInf,R_PosInf,NULL,NULL);
        double Mq=(np_plot_nn_zk(ker,a)-np_plot_nn_zk(ker,b))*v/h;
        Cq=(K*(1.0-u*v)+np_plot_nn_zk(ker,a)*v)/h/M-C*Mq/M;
        if(bw==BW_ADAP_NN) {
          f=K/h/M;fq=Ku/h/h/M;
        } else {
          double l=R_FINITE(lb)?(lb-q)/h:R_NegInf,r=R_FINITE(ub)?(ub-q)/h:R_PosInf;
          double T=np_cker_interval_mass(ker,l,r);
          if(!R_FINITE(T)||T<=0.0)error("np_plot_nn_unavailable: NN response query normalizer is unavailable");
          np_plot_nn_face(ker,l,-1.0-l*v);np_plot_nn_face(ker,r,-1.0-r*v);
          double Tq=(np_plot_nn_k(ker,l)+np_plot_nn_zk(ker,l)*v-
                     np_plot_nn_k(ker,r)-np_plot_nn_zk(ker,r)*v)/h;
          f=K/h/T;fq=Ku*(1.0-u*v)/h/h/T-f*(v/h+Tq/T);
        }
      }
      if(!R_FINITE(error_bound)||!R_FINITE(pdf_error)||!R_FINITE(pdfq_error)||
         error_bound<0.0||pdf_error<0.0||pdfq_error<0.0)
        error("np_plot_beta_numerical_unavailable: response discrepancy is not finite and nonnegative");
      if(!R_FINITE(C)||!R_FINITE(f)||!R_FINITE(Cq)||!R_FINITE(fq))
        error("np_plot_nn_unavailable: NN response jet is not representable");
      size_t ij=(size_t)j+(size_t)m*i;
      REAL(VECTOR_ELT(out,0))[ij]=C;REAL(VECTOR_ELT(out,1))[ij]=f;
      REAL(VECTOR_ELT(out,2))[ij]=Cq;REAL(VECTOR_ELT(out,3))[ij]=fq;
      REAL(VECTOR_ELT(out,4))[ij]=error_bound;
      REAL(VECTOR_ELT(out,8))[ij]=pdf_error;
      REAL(VECTOR_ELT(out,9))[ij]=pdfq_error;
    }
  }
  UNPROTECT(2);return out;
}
