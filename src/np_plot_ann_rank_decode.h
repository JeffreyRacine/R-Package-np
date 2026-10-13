/* New-request metadata only. Decode with the actual point owner's canonical
 * rounding/extended-NN option contract; never alter the user's stored counts. */
SEXP C_np_plot_ann_rank_decode(SEXP requested,SEXP nobs)
{
  if(TYPEOF(requested)!=REALSXP||XLENGTH(requested)<1||XLENGTH(requested)>INT_MAX||
     TYPEOF(nobs)!=INTSXP||XLENGTH(nobs)!=1||INTEGER(nobs)[0]==NA_INTEGER||INTEGER(nobs)[0]<2)
    error("invalid ANN rank decoder inputs");
  const int p=(int)XLENGTH(requested),n=INTEGER(nobs)[0];
  SEXP out=PROTECT(allocVector(VECSXP,3)),names=PROTECT(allocVector(STRSXP,3));
  SET_VECTOR_ELT(out,0,allocVector(INTSXP,p));SET_VECTOR_ELT(out,1,allocVector(REALSXP,p));SET_VECTOR_ELT(out,2,allocVector(LGLSXP,p));
  const char *keys[]={"lookup","distance.scale","extended"};
  for(int j=0;j<3;++j)SET_STRING_ELT(names,j,mkChar(keys[j]));setAttrib(out,R_NamesSymbol,names);
  for(int j=0;j<p;++j) {
    int k=0,extended=0;double scale=NA_REAL;
    if(np_nn_lookup_from_scale(n,1,REAL(requested)[j],&k,&scale,&extended))
      error("invalid canonical ANN rank under current options");
    INTEGER(VECTOR_ELT(out,0))[j]=k;REAL(VECTOR_ELT(out,1))[j]=scale;LOGICAL(VECTOR_ELT(out,2))[j]=extended;
  }
  UNPROTECT(2);return out;
}
