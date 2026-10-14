#include "np_plot_lp_precision_certificate.h"
/* Scratch request-only ordinary-regression beta moment-scale replay.
 * Include at end of jksum.c. No existing caller or numerical owner is edited.
 * The accepted B and P must come from the same ordinary LSQ response fit.
 */
typedef struct {
  SEXP args;
  int support_columns;
  const int *counts;
  const double *const *values;
  int previous_scale, scale_installed, derivative;
  int *terms, nterms, context_count;
  NPGLPBasisCtx *basis_context;
  NPLPSolveWorkspace solve;
  NPLPDesignSupport design;
} NPPlotBetaMomentCall;

static void np_plot_beta_moment_cleanup(void *data, Rboolean jump)
{
  NPPlotBetaMomentCall *c=(NPPlotBetaMomentCall *)data;
  (void)jump;
  if(c->scale_installed) int_LARGE_SF=c->previous_scale;
  np_lp_solve_workspace_clear(&c->solve);
  if(c->basis_context)for(int d=0;d<c->context_count;++d)
    np_glp_basis_ctx_free(c->basis_context+d);
  free(c->basis_context);free(c->terms);
}
static void np_plot_beta_moment_dims(SEXP x,int *n,int *p)
{
  SEXP d=getAttrib(x,R_DimSymbol);
  if(TYPEOF(x)!=REALSXP||TYPEOF(d)!=INTSXP||XLENGTH(d)!=2||
     INTEGER(d)[0]<1||INTEGER(d)[1]<0||
     XLENGTH(x)!=(R_xlen_t)INTEGER(d)[0]*INTEGER(d)[1])
    error("beta moment replay requires double matrices");
  *n=INTEGER(d)[0];*p=INTEGER(d)[1];
  for(R_xlen_t i=0;i<XLENGTH(x);++i)
    if(!R_FINITE(REAL(x)[i]))error("beta moment replay requires finite matrices");
}
static double **np_plot_beta_moment_cols(SEXP x,int n,int p)
{
  double **a=(double **)R_alloc(MAX(1,p),sizeof(double *));
  for(int d=0;d<p;++d)a[d]=REAL(x)+(size_t)d*n;
  return a;
}
static double **np_plot_beta_moment_workspace(int n,int p)
{
  double **a=(double **)R_alloc(MAX(1,p),sizeof(double *));
  for(int d=0;d<p;++d)a[d]=(double *)R_alloc(n,sizeof(double));
  return a;
}
static SEXP np_plot_beta_moment_body(void *data)
{
  NPPlotBetaMomentCall *c=(NPPlotBetaMomentCall *)data;
  static const char *labels[]={"train","eval","basis","response","bandwidth",
    "kernel","lower","upper","support","eval.train.index","max.bytes",
    "accepted.projection","degree","basis.code","bernstein","coordinate","order"};
  const int nargs=c->derivative?17:11;
  if(TYPEOF(c->args)!=VECSXP||XLENGTH(c->args)!=nargs)
    error("beta moment replay requires eleven named arguments");
  SEXP names=getAttrib(c->args,R_NamesSymbol);
  if(TYPEOF(names)!=STRSXP||XLENGTH(names)!=nargs)
    error("beta moment replay requires exact argument names");
  for(int z=0;z<nargs;++z)
    if(STRING_ELT(names,z)==NA_STRING||strcmp(CHAR(STRING_ELT(names,z)),labels[z]))
      error("beta moment replay argument order differs");
#define A(z) VECTOR_ELT(c->args,(z))
  SEXP train=A(0),eval=A(1),B=A(2),Y=A(3),bw=A(4),code=A(5);
  if(TYPEOF(train)!=VECSXP||XLENGTH(train)!=3||TYPEOF(eval)!=VECSXP||XLENGTH(eval)!=3||
     TYPEOF(code)!=INTSXP||XLENGTH(code)!=5)
    error("beta moment replay geometry or kernel metadata differs");
  int n=0,m=0,dims[3]={0,0,0};
  for(int g=0;g<3;++g){
    int nt,ne,pt,pe;
    np_plot_beta_moment_dims(VECTOR_ELT(train,g),&nt,&pt);
    np_plot_beta_moment_dims(VECTOR_ELT(eval,g),&ne,&pe);
    if(pt!=pe||(g&&(nt!=n||ne!=m)))error("beta moment replay geometry dimensions differ");
    n=nt;m=ne;dims[g]=pt;
  }
  const int nu=dims[0],no=dims[1],pc=dims[2];
  const double total=(double)nu+no+pc;
  int nb,w;np_plot_beta_moment_dims(B,&nb,&w);
  const int mode=INTEGER(code)[0],order=INTEGER(code)[1],
    ku=INTEGER(code)[2],ko=INTEGER(code)[3],compress=INTEGER(code)[4];
  if(n<2||pc<1||total>INT_MAX||nb!=n||w<1||w>INT_MAX-2||
     (mode!=BW_FIXED&&mode!=BW_GEN_NN&&mode!=BW_ADAP_NN)||
     ku<0||ku>=OP_UFUN_OFFSETS[1]||ko<0||ko>=OP_OFUN_OFFSETS[1]||
     (compress!=0&&compress!=1)||TYPEOF(bw)!=REALSXP||XLENGTH(bw)!=(R_xlen_t)total||
     TYPEOF(Y)!=REALSXP||XLENGTH(Y)!=n||
     TYPEOF(A(6))!=REALSXP||XLENGTH(A(6))!=pc||
     TYPEOF(A(7))!=REALSXP||XLENGTH(A(7))!=pc||
     TYPEOF(A(9))!=INTSXP||(XLENGTH(A(9))!=0&&XLENGTH(A(9))!=m)||
     TYPEOF(A(10))!=REALSXP||XLENGTH(A(10))!=1)
    error("beta moment replay invalid dimensions/options");
  np_continuous_kernel_descriptor descriptor;
  if(np_continuous_kernel_descriptor_init(NP_CKERNEL_FAMILY_BETA,
      NP_CKERNEL_COORDINATE_CODE,order,&descriptor)!=NP_CKERNEL_DESCRIPTOR_OK)
    error("beta moment replay invalid kernel order");
  const double forecast=64.0*(double)n*m+
    128.0*((double)n+m)*(total+w+3)+256.0*(double)w*(w+2)+16.0*(double)w*w+96.0*w+32.0;
  const double budget=REAL(A(10))[0];
  if(!R_FINITE(budget)||budget<=0||!R_FINITE(forecast)||forecast>budget||
     (double)n*m>INT_MAX||(double)w*(w+2)>INT_MAX)
    error("beta moment replay exceeds resource budget");
  for(int i=0;i<n;++i)
    if(!R_FINITE(REAL(Y)[i])||!R_FINITE(REAL(Y)[i]*REAL(Y)[i]))
      error("beta moment replay response square is not representable");
  for(int d=0;d<(int)total;++d)
    if(!R_FINITE(REAL(bw)[d])||(d<pc&&REAL(bw)[d]<=0))
      error("beta moment replay physical bandwidth is invalid");
  for(int d=0;d<pc;++d){
    const double lo=REAL(A(6))[d],hi=REAL(A(7))[d];
    if(!R_FINITE(lo)||!R_FINITE(hi)||!(lo<hi))error("beta moment replay needs finite support");
    for(int i=0;i<n;++i){double z=REAL(VECTOR_ELT(train,2))[(size_t)d*n+i];
      if(z<lo||z>hi)error("beta moment replay training row is outside support");}
    for(int i=0;i<m;++i){double z=REAL(VECTOR_ELT(eval,2))[(size_t)d*m+i];
      if(z<lo||z>hi)error("beta moment replay evaluation row is outside support");}
  }
  if(c->support_columns!=nu+no||(nu+no>0&&(c->counts==NULL||c->values==NULL)))
    error("beta moment replay categorical support differs");
  if(XLENGTH(A(9)))
    for(int j=0;j<m;++j){
      int i=INTEGER(A(9))[j];
      if(i<0||i>=n)error("beta moment replay occurrence map is invalid");
      for(int g=0;g<3;++g)for(int d=0;d<dims[g];++d)
        if(REAL(VECTOR_ELT(eval,g))[(size_t)d*m+j]!=
           REAL(VECTOR_ELT(train,g))[(size_t)d*n+i])
          error("beta moment replay occurrence map differs from training row");
    }
#ifdef MPI2
  if(!np_mpi_local_regression_active()&&iNum_Processors>1)
    error("beta moment replay requires explicit local regression scope");
#endif
  double **t[3],**e[3],**one[3];
  for(int g=0;g<3;++g){
    t[g]=np_plot_beta_moment_cols(VECTOR_ELT(train,g),n,dims[g]);
    e[g]=np_plot_beta_moment_cols(VECTOR_ELT(eval,g),m,dims[g]);
    one[g]=(double **)R_alloc(MAX(1,dims[g]),sizeof(double *));
  }
  double **basis=np_plot_beta_moment_cols(B,n,w);
  int which=0,derivative_order=0,bernstein=0,basis_mode=0;
  double *eval_basis=NULL,*eval_derivative=NULL,*coefficient=NULL,*adjoint=NULL,*adjoint_error=NULL;
  NPPlotLPCertificate certificate = {0};
  if(c->derivative){
    int np,wp;np_plot_beta_moment_dims(A(11),&np,&wp);
    if(np!=m||wp!=w||TYPEOF(A(12))!=INTSXP||XLENGTH(A(12))!=pc)
      error("beta derivative replay accepted projection or degree shape differs");
    for(int z=13;z<=16;++z)
      if(TYPEOF(A(z))!=INTSXP||XLENGTH(A(z))!=1||INTEGER(A(z))[0]==NA_INTEGER)
        error("beta derivative replay needs integer basis/direction metadata");
    basis_mode=INTEGER(A(13))[0];bernstein=INTEGER(A(14))[0];
    which=INTEGER(A(15))[0];derivative_order=INTEGER(A(16))[0];
    if((bernstein!=0&&bernstein!=1)||which<0||which>=pc||derivative_order<0||
       derivative_order>INTEGER(A(12))[which])
      error("beta derivative replay direction is unavailable");
    np_plot_validate_basis_request(pc,INTEGER(A(12)),basis_mode,w);
    if(!np_glp_build_terms(pc,INTEGER(A(12)),basis_mode,&c->terms,&c->nterms)||c->nterms!=w)
      error("beta derivative replay canonical term width differs");
    if(bernstein){
      c->basis_context=(NPGLPBasisCtx *)calloc(pc,sizeof(NPGLPBasisCtx));
      if(!c->basis_context)error("beta derivative replay basis allocation failed");
      c->context_count=pc;
      for(int d=0;d<pc;++d){
        double lo=t[2][d][0],hi=lo;
        for(int i=1;i<n;++i){if(t[2][d][i]<lo)lo=t[2][d][i];if(t[2][d][i]>hi)hi=t[2][d][i];}
        if(!np_glp_basis_ctx_init(c->basis_context+d,INTEGER(A(12))[d],lo,hi,basis_mode==1))
          error("beta derivative replay canonical basis initialization failed");
      }
    }
    double **rebuilt=np_plot_beta_moment_workspace(n,w);
    if(bernstein)np_glp_fill_basis_train(pc,c->terms,w,t[2],n,c->basis_context,rebuilt);
    else np_glp_fill_basis_raw_train(pc,c->terms,w,t[2],n,rebuilt);
    for(int k=0;k<w;++k)for(int i=0;i<n;++i)
      if(rebuilt[k][i]!=basis[k][i])error("beta derivative replay basis differs from accepted fit");
    np_lp_design_support_prepare(&c->design,basis,n,w);
    if(!np_lp_solve_workspace_reserve(&c->solve,w,1))
      error("beta derivative replay solve allocation failed");
    eval_basis=(double *)R_alloc(w,sizeof(double));
    eval_derivative=(double *)R_alloc(w,sizeof(double));
    coefficient=(double *)R_alloc(w,sizeof(double));
    adjoint=(double *)R_alloc(w,sizeof(double));adjoint_error=(double *)R_alloc(w,sizeof(double));
    np_plot_lp_certificate_init(&certificate,w);
  }
  double *y2=(double *)R_alloc(n,sizeof(double));
  double **outer=(double **)R_alloc(w+2,sizeof(double *));
  for(int i=0;i<n;++i)y2[i]=REAL(Y)[i]*REAL(Y)[i];
  outer[0]=y2;outer[1]=REAL(Y);for(int k=0;k<w;++k)outer[k+2]=basis[k];
  double *moments=(double *)R_alloc((size_t)w*(w+2),sizeof(double));
  int *op=(int *)R_alloc((size_t)total,sizeof(int));
  for(int k=0;k<(int)total;++k)op[k]=OP_NORMAL;
  int *kernel_u=(int *)R_alloc(MAX(1,nu),sizeof(int));
  int *kernel_o=(int *)R_alloc(MAX(1,no),sizeof(int));
  for(int k=0;k<nu;++k)kernel_u[k]=ku;
  for(int k=0;k<no;++k)kernel_o[k]=ko;
  int *counts=(int *)R_alloc(MAX(1,nu+no),sizeof(int));
  double **values=(double **)R_alloc(MAX(1,nu+no),sizeof(double *));
  for(int k=0;k<nu+no;++k){
    counts[k]=c->counts[k];values[k]=(double *)c->values[k];
  }
  int nr=mode==BW_FIXED?1:(mode==BW_GEN_NN?m:n);
  double **h=np_plot_beta_moment_workspace(nr,pc);
  double **he=(double **)R_alloc(pc,sizeof(double *));
  double *lambda=(double *)R_alloc(MAX(1,nu+no),sizeof(double));
  NPNNGeometryContext geometry={.mode=XLENGTH(A(9))?NP_NN_QUERY_TRAINING_MAP:NP_NN_QUERY_EXTERNAL,
    .eval_to_train=XLENGTH(A(9))?INTEGER(A(9)):NULL,.adaptive_successor=NULL};
  NPNNGeometryStatus status=NP_NN_GEOMETRY_OK;
  c->previous_scale=int_LARGE_SF;c->scale_installed=1;int_LARGE_SF=SF_ARB;
  if(np_beta_continuous_bandwidth_prepare_canonical_ctx(mode,n,m,nu,no,pc,
      t[2],e[2],REAL(bw),h,NULL,lambda,NULL,&geometry,&status))
    error("beta moment replay canonical bandwidth preparation failed (%d)",status);
  NPContinuousKernelRoute route={0};route.segment_count=1;
  route.segment[0].descriptor=descriptor;route.segment[0].coordinate_offset=0;
  route.segment[0].coordinate_count=pc;route.segment[0].lower=REAL(A(6));route.segment[0].upper=REAL(A(7));
  NPContinuousKernelDerivativeDiagnostics diag={.bad_coordinate=-1,.bad_observation=-1,
    .undefined_count=0,.beta_status=NP_BETA_OK};
  int nprotect=1;
  SEXP result=PROTECT(allocMatrix(REALSXP,n,m));
  SEXP D=R_NilValue,P=R_NilValue,point=R_NilValue,gradient=R_NilValue,ridge=R_NilValue,map_error=R_NilValue,relative_error=R_NilValue,defect=R_NilValue,certified=R_NilValue;
  if(c->derivative){
    D=PROTECT(allocMatrix(REALSXP,n,m));P=PROTECT(allocMatrix(REALSXP,m,w));
    point=PROTECT(allocVector(REALSXP,m));gradient=PROTECT(allocVector(REALSXP,m));
    ridge=PROTECT(allocVector(REALSXP,m));map_error=PROTECT(allocMatrix(REALSXP,n,m));
    relative_error=PROTECT(allocVector(REALSXP,m));defect=PROTECT(allocVector(REALSXP,m));
    certified=PROTECT(allocVector(LGLSXP,m));nprotect+=9;
  }
  for(int j=0;j<m;++j){
    R_CheckUserInterrupt();
    for(int g=0;g<3;++g)for(int d=0;d<dims[g];++d)one[g][d]=e[g][d]+j;
    for(int d=0;d<pc;++d)he[d]=mode==BW_GEN_NN?h[d]+j:h[d];
    NPLPSupportRankContext rank={0,w,c->derivative?&c->design:NULL};
    if(np_beta_regression_lp_moment_row_canonical(mode,n,nu,no,pc,w,w+2,
       op,kernel_u,kernel_o,t[0],t[1],t[2],one[0],one[1],one[2],
       outer,basis,REAL(bw),h,he,(nu+no)>0?lambda:NULL,counts,values,compress,
       moments,NULL,NULL,REAL(result)+(size_t)n*j,&route,&diag,c->derivative?&rank:NULL))
      error("beta moment replay canonical moment row failed (%d)",diag.beta_status);
    for(int i=0;i<n;++i)
      if(!R_FINITE(REAL(result)[(size_t)n*j+i]))error("beta moment replay produced nonfinite weights");
    if(c->derivative){
      for(int i=0;i<w;++i){
        c->solve.rhs_source[i]=moments[i*(w+2)+1];
        for(int k=0;k<w;++k)c->solve.gram_source[i+(size_t)k*w]=moments[i*(w+2)+2+k];
      }
      const double anchor=c->solve.gram_source[0];
      NPLPSolvePolicyDiagnostics policy={0,0.0};
      if(np_lp_solve_workspace_solve_response_ranked(&c->solve,w,1,1.0/n,rank.count,&policy)!=NP_LP_SOLVE_POLICY_OK)
        error("beta derivative replay response owner did not accept row");
      for(int i=0;i<w;++i)coefficient[i]=c->solve.rhs_work[i];
      if(bernstein){
        np_glp_fill_basis_eval(pc,c->terms,w,e[2],j,c->basis_context,eval_basis);
        if(derivative_order>0)np_glp_fill_basis_eval_deriv(which,derivative_order,pc,c->terms,w,e[2],j,c->basis_context,eval_derivative);
        else memcpy(eval_derivative,eval_basis,w*sizeof(double));
      }else{
        np_glp_fill_basis_eval_raw(pc,c->terms,w,e[2],j,eval_basis);
        if(derivative_order>0)np_glp_fill_basis_eval_deriv_raw(which,derivative_order,pc,c->terms,w,e[2],j,eval_derivative);
        else memcpy(eval_derivative,eval_basis,w*sizeof(double));
      }
      REAL(point)[j]=0.0;REAL(gradient)[j]=0.0;
      for(int i=0;i<w;++i){
        REAL(point)[j]+=eval_basis[i]*coefficient[i];
        REAL(gradient)[j]+=eval_derivative[i]*coefficient[i];
        c->solve.rhs_source[i]=eval_basis[i];
      }
      if(np_lp_solve_workspace_solve_adjoint_factored(&c->solve,w,1,anchor,&policy)!=NP_LP_SOLVE_POLICY_OK)
        error("beta derivative replay level adjoint failed");
      for(int i=0;i<w;++i){
        const double v=c->solve.rhs_work[i],accepted=REAL(A(11))[(size_t)i*m+j];
        if(!R_FINITE(v)||fabs(v-accepted)>2e-11*(1.0+fabs(accepted)))
          error("beta derivative replay level adjoint differs from accepted fit");
        REAL(P)[(size_t)i*m+j]=v;
        c->solve.rhs_source[i]=eval_derivative[i];
      }
      if(np_lp_solve_workspace_solve_adjoint_factored(&c->solve,w,1,anchor,&policy)!=NP_LP_SOLVE_POLICY_OK)
        error("beta derivative replay derivative adjoint failed");
      for(int k=0;k<w;++k)adjoint[k]=derivative_order==0 ?
        REAL(A(11))[j+(size_t)k*m] : c->solve.rhs_work[k];
      for(int i=0;i<n;++i){
        double v=0.0;for(int k=0;k<w;++k)v+=basis[k][i]*adjoint[k];
        v*=REAL(result)[(size_t)n*j+i];
        if(!R_FINITE(v))error("beta derivative replay derivative map is not representable");
        REAL(D)[(size_t)n*j+i]=v;
      }
      REAL(ridge)[j]=policy.ridge_total;
      double defect_value=R_PosInf;
      int precision=np_plot_lp_certificate_adjoint(&certificate,&c->solve,eval_derivative,
        adjoint,anchor,&policy,adjoint_error,&defect_value);
      REAL(defect)[j]=defect_value;
      double map_scale=0.0,error_scale=0.0;
      for(int i=0;i<n;++i){
        const size_t at=(size_t)n*j+i;
        const double err=precision?np_plot_lp_certificate_map_error(w,basis,i,
          REAL(result)[at],adjoint,adjoint_error,REAL(D)[at]):R_PosInf;
        REAL(map_error)[at]=err;
        if(!R_FINITE(err))precision=0;
        if(fabs(REAL(D)[at])>map_scale)map_scale=fabs(REAL(D)[at]);
        if(err>error_scale)error_scale=err;
      }
      double ratio=R_PosInf;
      if(precision&&error_scale==0.0)ratio=0.0;
      else if(precision&&map_scale>0.0){
        double ss_map=0.0,ss_error=0.0;
        for(int i=0;i<n;++i){const size_t at=(size_t)n*j+i;
          const double a=nextafter(fabs(REAL(D)[at])/map_scale,0.0);
          const double b=REAL(map_error)[at]==0.0?0.0:nextafter(REAL(map_error)[at]/error_scale,INFINITY);
          ss_map=fmax(0.0,nextafter(ss_map+fmax(0.0,nextafter(a*a,0.0)),0.0));
          ss_error=np_plot_up_add(ss_error,np_plot_up_mul(b,b));}
        if(ss_map>0.0)ratio=np_plot_up_mul(nextafter(error_scale/map_scale,INFINITY),
          nextafter(sqrt(nextafter(ss_error/ss_map,INFINITY)),INFINITY));
      }
      REAL(relative_error)[j]=ratio;
      LOGICAL(certified)[j]=precision&&R_FINITE(ratio)&&ratio<=2e-11;
    }
  }
  /* Root the new value before symbol installation, regardless of the C
   * compiler's argument evaluation order. */
  SEXP forecast_value=PROTECT(ScalarReal(forecast));nprotect++;
  setAttrib(result,install("forecast.bytes"),forecast_value);
  if(!c->derivative){UNPROTECT(nprotect);return result;}
  SEXP out=PROTECT(allocVector(VECSXP,11)),outnames=PROTECT(allocVector(STRSXP,11));nprotect+=2;
  SEXP requested_order=PROTECT(ScalarInteger(derivative_order));nprotect++;
  SEXP value[11]={result,D,P,point,gradient,ridge,map_error,relative_error,defect,certified,requested_order};
  const char *keys[]={"kernel.weights","derivative.weights","level.projection","mean","gradient","ridge","operator.error","operator.relative.error","inverse.defect","certified","requested.order"};
  for(int i=0;i<11;++i){SET_VECTOR_ELT(out,i,value[i]);SET_STRING_ELT(outnames,i,mkChar(keys[i]));}
  setAttrib(out,R_NamesSymbol,outnames);UNPROTECT(nprotect);return out;
#undef A
}
SEXP np_plot_beta_moment_weights_native(SEXP args,int support_columns,
  const int *counts,const double *const *values)
{
  NPPlotBetaMomentCall c={.args=args,.support_columns=support_columns,
    .counts=counts,.values=values,.scale_installed=0,.derivative=TYPEOF(args)==VECSXP&&XLENGTH(args)==17};
  np_lp_solve_workspace_init(&c.solve);
  return R_UnwindProtect(np_plot_beta_moment_body,&c,np_plot_beta_moment_cleanup,&c,NULL);
}
