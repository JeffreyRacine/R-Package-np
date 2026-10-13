/* Request-only LSQ inference bridge. Include in np.c after the canonical
 * descriptor parser and np_regression declaration. Existing wrappers unchanged. */
static int np_plot_projection_scalar_int(SEXP x, const char *what)
{
  if(TYPEOF(x) != INTSXP || XLENGTH(x) != 1 || INTEGER(x)[0] == NA_INTEGER)
    error("plot gradient projection: invalid %s", what);
  return INTEGER(x)[0];
}
static void np_plot_projection_real_length(SEXP x, R_xlen_t length,
                                           const char *what, int finite)
{
  if(TYPEOF(x) != REALSXP || XLENGTH(x) != length)
    error("plot gradient projection: invalid %s length/type", what);
  if(finite)
    for(R_xlen_t i=0;i<length;++i)
      if(!R_FINITE(REAL(x)[i]))
        error("plot gradient projection: nonfinite %s", what);
}
SEXP C_np_plot_regression_gradient_projection(SEXP args)
{
  if(TYPEOF(args)!=VECSXP || XLENGTH(args)!=24)
    error("plot gradient projection requires 24 typed arguments");
  SEXP v[24];
  for(int i=0;i<24;++i)v[i]=VECTOR_ELT(args,i);
  if(TYPEOF(v[14])!=INTSXP || XLENGTH(v[14])!=REG_CKORDERI+1)
    error("plot gradient projection: invalid option vector");
  const int *op=INTEGER(v[14]);
  for(int i=0;i<=REG_CKORDERI;++i)
    if(op[i]==NA_INTEGER)error("plot gradient projection: missing option");
  const int n=op[REG_TNOBSI],m=op[REG_ENOBSI];
  const int nu=op[REG_NUNOI],no=op[REG_NORDI],pc=op[REG_NCONI];
  const int en=np_plot_projection_scalar_int(v[19],"evaluation count");
  const int nc=np_plot_projection_scalar_int(v[20],"predictor count");
  int bern=np_plot_projection_scalar_int(v[17],"Bernstein flag");
  int basis=np_plot_projection_scalar_int(v[18],"basis code");
  if(n<2 || m<1 || pc<1 || nu<0 || no<0 ||
     (double)nu+no+pc>INT_MAX || nc!=nu+no+pc || en!=m ||
     op[REG_LL]!=NP_LP_ENGINE_GENERAL || op[REG_GRAD]!=1 ||
     op[REG_LSFI]!=SF_ARB || op[REG_EY]!=1 ||
     (op[REG_TISEI]!=0 && op[REG_TISEI]!=1) ||
     (op[REG_TISEI] && n!=m) || (bern!=0 && bern!=1) ||
     (op[REG_CATCOMPI]!=0 && op[REG_CATCOMPI]!=1))
    error("plot gradient projection: inconsistent dimensions or request");
  if(TYPEOF(v[21])!=INTSXP || XLENGTH(v[21])!=3 ||
     INTEGER(v[21])[0]!=NP_REGRESSION_OUTPUT_MEAN_AND_GRADIENT ||
     INTEGER(v[21])[1]!=0 || INTEGER(v[21])[2]<1 ||
     INTEGER(v[21])[2]==NA_INTEGER ||
     getAttrib(v[21],install(".np.gradient.coordinate"))!=R_NilValue)
    error("plot gradient projection: expected gradient-only compact projection request");
  const int width=INTEGER(v[21])[2];
  const double resident=(double)width*(n+(double)m)+
    2.0*(n+(double)m)*nc+2.0*n+2.0*m+(double)m*nc+256.0*width*width;
  if(!R_FINITE(resident) || resident>8.0*1024*1024)
    error("plot gradient projection: allocation budget exceeded");
  if(TYPEOF(v[15])!=INTSXP || XLENGTH(v[15])!=pc ||
     TYPEOF(v[16])!=INTSXP || XLENGTH(v[16])!=pc)
    error("plot gradient projection: invalid degree or derivative metadata");
  int any_positive=0,any_available=0;
  for(int i=0;i<pc;++i){
    const int d=INTEGER(v[15])[i],s=INTEGER(v[16])[i];
    if(d<0 || s<1 || d==NA_INTEGER || s==NA_INTEGER)
      error("plot gradient projection: invalid degree/order");
    any_positive|=d>0;any_available|=s<=d;
  }
  if(basis<0 || basis>2 ||
     (op[REG_DOTREEI]!=NP_TREE_FALSE && op[REG_DOTREEI]!=NP_TREE_TRUE) ||
     (op[REG_MINIOI]!=IO_MIN_FALSE && op[REG_MINIOI]!=IO_MIN_TRUE) ||
     op[REG_UKRNEVI]<0 || op[REG_UKRNEVI]>=OP_UFUN_OFFSETS[1] ||
     op[REG_OKRNEVI]<0 || op[REG_OKRNEVI]>=OP_OFUN_OFFSETS[1])
    error("plot gradient projection: invalid basis/tree/I-O/category-kernel options");
  np_plot_validate_basis_request(pc,INTEGER(v[15]),basis,width);
  if(!any_positive || !any_available)
    error("plot gradient projection: no available positive-degree gradient");
  np_plot_projection_real_length(v[0],(R_xlen_t)n*nu,"training unordered",1);
  np_plot_projection_real_length(v[1],(R_xlen_t)n*no,"training ordered",1);
  np_plot_projection_real_length(v[2],(R_xlen_t)n*pc,"training continuous",1);
  np_plot_projection_real_length(v[3],n,"response",1);
  const int ne=op[REG_TISEI]?0:m;
  np_plot_projection_real_length(v[4],(R_xlen_t)ne*nu,"evaluation unordered",1);
  np_plot_projection_real_length(v[5],(R_xlen_t)ne*no,"evaluation ordered",1);
  np_plot_projection_real_length(v[6],(R_xlen_t)ne*pc,"evaluation continuous",1);
  np_plot_projection_real_length(v[7],0,"unused evaluation response",1);
  np_plot_projection_real_length(v[8],nc,"physical bandwidth",1);
  if(op[REG_MLEVI]<0 || ((nu+no)>0 && op[REG_MLEVI]<1))
    error("plot gradient projection: invalid declared level shape");
  np_plot_projection_real_length(v[9],(R_xlen_t)op[REG_MLEVI]*(nu+no),"declared levels",0);
  for(int i=10;i<=12;++i)
    np_plot_projection_real_length(v[i],1,"scale/padding metadata",0);
  np_plot_projection_real_length(v[13],pc,"continuous scale",1);
  np_plot_projection_real_length(v[22],pc,"lower bounds",0);
  np_plot_projection_real_length(v[23],pc,"upper bounds",0);
  for(int i=0;i<pc;++i)
    if(ISNAN(REAL(v[22])[i]) || ISNAN(REAL(v[23])[i]) ||
       !(REAL(v[22])[i]<REAL(v[23])[i]))
      error("plot gradient projection: invalid support bounds");
  if(op[REG_BWI]!=BW_FIXED && op[REG_BWI]!=BW_GEN_NN && op[REG_BWI]!=BW_ADAP_NN)
    error("plot gradient projection: invalid bandwidth mode");
  np_continuous_kernel_descriptor descriptor=np_continuous_kernel_descriptor_or_error(
    op[REG_CKFAMILYI],op[REG_CKRNEVI],op[REG_CKORDERI],"plot gradient projection");
#ifdef MPI2
  if(!np_mpi_local_regression_active() && iNum_Processors>1)
    error("plot gradient projection requires an explicit local regression scope");
#endif
  /* The canonical owner may reorder or otherwise mutate its input arrays. */
  int protected_count=0;
  for(int i=0;i<17;++i){PROTECT(v[i]=duplicate(v[i]));++protected_count;}
  PROTECT(v[22]=duplicate(v[22]));++protected_count;
  PROTECT(v[23]=duplicate(v[23]));++protected_count;
  SEXP mean=PROTECT(allocVector(REALSXP,m));++protected_count;
  SEXP grad=PROTECT(allocVector(REALSXP,(R_xlen_t)m*nc));++protected_count;
  SEXP B=PROTECT(allocMatrix(REALSXP,n,width));++protected_count;
  SEXP P=PROTECT(allocMatrix(REALSXP,m,width));++protected_count;
  double extra[6]={0};int variance_unavailable=0;int *gradient_zero=NULL;
  NPRegressionFitProjection projection={.ntrain=n,.neval=m,.width=width,
    .basis=REAL(B),.projection=REAL(P)};
  NPContinuousKernelRoute beta_route;
  NPContinuousKernelDerivativeDiagnostics diagnostics;
  const NPContinuousKernelRoute *route=NULL;
  NPContinuousKernelDerivativeDiagnostics *diag=NULL;
  if(descriptor.family==NP_CKERNEL_FAMILY_BETA){
    beta_route.segment_count=1;
    beta_route.segment[0].descriptor=descriptor;
    beta_route.segment[0].coordinate_offset=0;
    beta_route.segment[0].coordinate_count=pc;
    beta_route.segment[0].lower=REAL(v[22]);beta_route.segment[0].upper=REAL(v[23]);
    diagnostics.bad_coordinate=-1;diagnostics.bad_observation=-1;
    diagnostics.undefined_count=0;diagnostics.beta_status=NP_BETA_OK;
    route=&beta_route;diag=&diagnostics;
  }
  np_regression(REAL(v[0]),REAL(v[1]),REAL(v[2]),REAL(v[3]),
    REAL(v[4]),REAL(v[5]),REAL(v[6]),REAL(v[7]),REAL(v[8]),REAL(v[9]),
    REAL(v[10]),REAL(v[11]),REAL(v[12]),REAL(v[13]),INTEGER(v[14]),
    INTEGER(v[15]),INTEGER(v[16]),&bern,&basis,REAL(mean),NULL,REAL(grad),NULL,
    extra,REAL(v[22]),REAL(v[23]),route,diag,op[REG_CATCOMPI],NULL,NULL,
    &variance_unavailable,&gradient_zero,&projection);
  if(projection.emitted!=m)
    error("plot gradient projection did not visit every evaluation row");
  SEXP projection_result=PROTECT(allocVector(VECSXP,3));++protected_count;
  SET_VECTOR_ELT(projection_result,0,B);SET_VECTOR_ELT(projection_result,1,P);
  SET_VECTOR_ELT(projection_result,2,ScalarLogical(projection.constant));
  SEXP out=PROTECT(allocVector(VECSXP,3));++protected_count;
  SET_VECTOR_ELT(out,0,mean);SET_VECTOR_ELT(out,1,grad);SET_VECTOR_ELT(out,2,projection_result);
  SEXP names=PROTECT(allocVector(STRSXP,3));++protected_count;
  SET_STRING_ELT(names,0,mkChar("mean"));SET_STRING_ELT(names,1,mkChar("gradient"));
  SET_STRING_ELT(names,2,mkChar("fit.projection"));setAttrib(out,R_NamesSymbol,names);
  UNPROTECT(protected_count);return out;
}
