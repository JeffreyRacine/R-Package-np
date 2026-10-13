/* Include in np.c after existing categorical support parsers. */
extern SEXP np_plot_beta_moment_weights_native(SEXP,int,const int *,const double *const *);
SEXP C_np_plot_beta_moment_weights(SEXP args)
{
  if(TYPEOF(args)!=VECSXP||(XLENGTH(args)!=11&&XLENGTH(args)!=17))
    error("beta moment replay requires eleven or seventeen arguments");
  SEXP train=VECTOR_ELT(args,0),support=VECTOR_ELT(args,8);
  if(TYPEOF(train)!=VECSXP||XLENGTH(train)!=3)
    error("beta moment replay requires three training groups");
  int n=0,dims[4]={0,0,0,0};
  SEXP data[4]={R_NilValue,R_NilValue,VECTOR_ELT(train,0),VECTOR_ELT(train,1)};
  for(int g=0;g<3;++g){
    int nr=0,nc=0;SEXP x=VECTOR_ELT(train,g);
    np_matrix_dims_or_zero(x,&nr,&nc);
    if(TYPEOF(x)!=REALSXP||nr<1||nc<0||XLENGTH(x)!=(R_xlen_t)nr*nc||(g&&nr!=n))
      error("beta moment replay training shape differs");
    n=nr;if(g<2)dims[g+2]=nc;
  }
  const NPCategoricalSupport s=np_categorical_support_parse(support,n,dims,data);
  return np_plot_beta_moment_weights_native(args,s.columns,s.counts,s.values);
}
