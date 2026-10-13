/* New-request ingress only. Bound canonical enumeration before allocating
 * terms, then ask the unchanged canonical builder to certify exact width. */
void np_plot_validate_basis_request(int p, const int *degree, int basis, int width)
{
  if(p<1 || p>1024 || degree==NULL || width<1 || basis<0 || basis>2 ||
     4.0*p*width>16.0*1024*1024)
    error("invalid plot projection basis request");
  for(int j=0;j<p;++j)if(degree[j]<0 || degree[j]>100)
    error("invalid plot projection degree");
  double count=1.0, enumeration=1.0;
  if(basis==1) {
    int max_degree=0;
    for(int j=0;j<p;++j)max_degree=MAX(max_degree,degree[j]);
    double previous[101]={1.0},next[101];
    for(int j=0;j<p;++j) {
      memset(next,0,sizeof(next));
      for(int total=0;total<=max_degree;++total)
        for(int k=0;k<=MIN(degree[j],total);++k)
          next[total]=MIN((double)width+1,next[total]+previous[total-k]);
      memcpy(previous,next,sizeof(previous));
    }
    count=0.0;for(int k=0;k<=max_degree;++k)count+=previous[k];
  } else for(int j=0;j<p;++j) {
    if(basis==0)count+=degree[j];else count*=1.0+degree[j];
    enumeration*=1.0+degree[j];
    if(count>width || enumeration>1e7)
      error("plot projection basis exceeds its request budget");
  }
  if(count!=width)error("plot projection basis width does not match its degree contract");
  int *terms=NULL,actual=0;
  const int ok=np_glp_build_terms(p,degree,basis,&terms,&actual);
  free(terms);
  if(!ok || actual!=width)error("canonical plot projection basis width mismatch");
}

/* New plot-inference request only. No existing fit calls this replay. */
typedef struct {
  SEXP basis, x, ex, degree, permutation, accepted, output;
  int n, p, m, width, bernstein, basis_mode, coordinate, order;
  int *terms;
  NPGLPBasisCtx *contexts;
  NPLPFullRowWorkspace inverse;
  double **eval_columns;
  double *level, *direction, *adjoint;
} NPPlotAllLargeDerivative;

static void np_plot_alllarge_derivative_cleanup(void *data, Rboolean jump)
{
  NPPlotAllLargeDerivative *c = (NPPlotAllLargeDerivative *)data;
  (void)jump;
  if(c->contexts != NULL) {
    for(int j=0;j<c->p;++j) np_glp_basis_ctx_free(&c->contexts[j]);
    free(c->contexts);
  }
  free(c->terms);
  np_lp_full_row_workspace_clear(&c->inverse);
}

static SEXP np_plot_alllarge_derivative_execute(void *data)
{
  NPPlotAllLargeDerivative *c=(NPPlotAllLargeDerivative *)data;
  int width=0;
  if(!np_glp_build_terms(c->p,INTEGER(c->degree),c->basis_mode,&c->terms,&width) || width!=c->width)
    error("plot derivative basis contract does not match the accepted owner");
  if(!np_lp_full_row_workspace_reserve(&c->inverse,width,1))
    error("plot derivative inverse workspace unavailable");
  if(c->bernstein) {
    c->contexts=(NPGLPBasisCtx *)calloc((size_t)c->p,sizeof(NPGLPBasisCtx));
    if(c->contexts==NULL) error("plot derivative basis workspace unavailable");
    for(int j=0;j<c->p;++j) {
      double lo=REAL(c->x)[c->n*j],hi=lo;
      for(int i=1;i<c->n;++i) {double v=REAL(c->x)[i+c->n*j];lo=MIN(lo,v);hi=MAX(hi,v);}
      if(!np_glp_basis_ctx_init(&c->contexts[j],INTEGER(c->degree)[j],lo,hi,c->basis_mode==1))
        error("plot derivative basis initialization failed");
    }
  }
  for(int a=0;a<width;++a)for(int b=0;b<width;++b)c->inverse.matrix_copy[a+b*width]=0.0;
  for(int k=0;k<c->n;++k) {
    int i=INTEGER(c->permutation)[k];
    for(int a=0;a<width;++a) {
      double za=REAL(c->basis)[i+c->n*a];
      for(int b=a;b<width;++b) {
        double v=za*REAL(c->basis)[i+c->n*b];
        c->inverse.matrix_copy[a+b*width]+=v;
        if(a!=b)c->inverse.matrix_copy[b+a*width]+=v;
      }
    }
  }
  double initial=c->inverse.matrix_copy[0];
  if(!np_lp_full_row_workspace_invert_retryable(&c->inverse,width,1.0/c->n,64))
    error("plot derivative accepted inverse replay failed");
  REAL(VECTOR_ELT(c->output,2))[0]=c->inverse.matrix_copy[0]-initial;
  if(!np_lp_full_row_workspace_pack_inverse_rows(&c->inverse,width))
    error("plot derivative inverse packing failed");
  for(int j=0;j<c->m;++j) {
    if(c->bernstein) {
      np_glp_fill_basis_eval(c->p,c->terms,width,c->eval_columns,j,c->contexts,c->level);
      np_glp_fill_basis_eval_deriv(c->coordinate,c->order,c->p,c->terms,width,c->eval_columns,j,c->contexts,c->direction);
    } else {
      np_glp_fill_basis_eval_raw(c->p,c->terms,width,c->eval_columns,j,c->level);
      np_glp_fill_basis_eval_deriv_raw(c->coordinate,c->order,c->p,c->terms,width,c->eval_columns,j,c->direction);
    }
    for(int k=0;k<width;++k) {
      double u=0.0,v=0.0;
      for(int l=0;l<width;++l) {
        u+=c->inverse.matrix_copy[l*width+k]*c->direction[l];
        v+=c->inverse.matrix_copy[l*width+k]*c->level[l];
      }
      double expected=REAL(c->accepted)[j+c->m*k];
      if(!R_FINITE(expected)||!R_FINITE(u)||!R_FINITE(v)||fabs(v-expected)>2e-11*(1.0+fabs(expected)))
        error("plot derivative replay does not reproduce the accepted projection");
      c->adjoint[k]=u;REAL(VECTOR_ELT(c->output,1))[j+c->m*k]=v;
    }
    for(int i=0;i<c->n;++i) {
      double h=0.0;
      for(int k=0;k<width;++k)h+=REAL(c->basis)[i+c->n*k]*c->adjoint[k];
      if(!R_FINITE(h))error("plot derivative response map is non-finite");
      REAL(VECTOR_ELT(c->output,0))[j+c->m*i]=h;
    }
  }
  return c->output;
}

SEXP C_np_plot_alllarge_derivative(SEXP basis, SEXP x, SEXP ex, SEXP degree,
  SEXP bernstein, SEXP basis_mode, SEXP coordinate, SEXP order, SEXP permutation, SEXP accepted)
{
  SEXP bd=getAttrib(basis,R_DimSymbol),xd=getAttrib(x,R_DimSymbol),ed=getAttrib(ex,R_DimSymbol),ad=getAttrib(accepted,R_DimSymbol);
  if(TYPEOF(basis)!=REALSXP||TYPEOF(x)!=REALSXP||TYPEOF(ex)!=REALSXP||TYPEOF(accepted)!=REALSXP||
     TYPEOF(bd)!=INTSXP||TYPEOF(xd)!=INTSXP||TYPEOF(ed)!=INTSXP||TYPEOF(ad)!=INTSXP||
     XLENGTH(bd)!=2||XLENGTH(xd)!=2||XLENGTH(ed)!=2||XLENGTH(ad)!=2||TYPEOF(degree)!=INTSXP||TYPEOF(permutation)!=INTSXP)
    error("invalid plot derivative matrix contract");
  if(XLENGTH(coordinate)!=1||XLENGTH(order)!=1||XLENGTH(bernstein)!=1||XLENGTH(basis_mode)!=1||
     asInteger(coordinate)==NA_INTEGER||asInteger(order)==NA_INTEGER||asInteger(basis_mode)==NA_INTEGER)
    error("invalid plot derivative scalar contract");
  NPPlotAllLargeDerivative c={.basis=basis,.x=x,.ex=ex,.degree=degree,.permutation=permutation,.accepted=accepted,
    .n=INTEGER(bd)[0],.width=INTEGER(bd)[1],.p=INTEGER(xd)[1],.m=INTEGER(ed)[0],
    .bernstein=asLogical(bernstein),.basis_mode=asInteger(basis_mode),.coordinate=asInteger(coordinate)-1,.order=asInteger(order)};
  if(c.n<1||c.p<1||c.m<1||c.width<1||INTEGER(xd)[0]!=c.n||INTEGER(ed)[1]!=c.p||INTEGER(ad)[0]!=c.m||INTEGER(ad)[1]!=c.width||
     LENGTH(degree)!=c.p||LENGTH(permutation)!=c.n||c.coordinate<0||c.coordinate>=c.p||c.order<1||c.order>100||
     c.bernstein==NA_LOGICAL||c.basis_mode<0||c.basis_mode>2)
    error("invalid plot derivative dimensions or direction");
  if(XLENGTH(basis)!=(R_xlen_t)c.n*c.width||XLENGTH(x)!=(R_xlen_t)c.n*c.p||
     XLENGTH(ex)!=(R_xlen_t)c.m*c.p||XLENGTH(accepted)!=(R_xlen_t)c.m*c.width)
    error("invalid plot derivative matrix extent");
  double bytes=8.0*(20.0*c.width*c.width+2.0*c.m*c.width+(double)c.m*c.n+8.0*c.width+4.0*c.p*c.width)+4.0*c.n;
  if(!R_FINITE(bytes)||bytes>64.0*1024*1024)error("plot derivative request exceeds its workspace budget");
  np_plot_validate_basis_request(c.p,INTEGER(degree),c.basis_mode,c.width);
  for(R_xlen_t i=0;i<XLENGTH(basis);++i)if(!R_FINITE(REAL(basis)[i]))error("non-finite plot derivative basis");
  for(R_xlen_t i=0;i<XLENGTH(x);++i)if(!R_FINITE(REAL(x)[i]))error("non-finite plot derivative training geometry");
  for(R_xlen_t i=0;i<XLENGTH(ex);++i)if(!R_FINITE(REAL(ex)[i]))error("non-finite plot derivative evaluation geometry");
  int *seen=(int *)R_alloc((size_t)c.n,sizeof(int));memset(seen,0,(size_t)c.n*sizeof(int));
  for(int i=0;i<c.n;++i){int k=INTEGER(permutation)[i];if(k<0||k>=c.n||seen[k])error("invalid plot derivative donor permutation");seen[k]=1;}
  c.eval_columns=(double **)R_alloc((size_t)c.p,sizeof(double *));
  for(int j=0;j<c.p;++j)c.eval_columns[j]=REAL(ex)+(size_t)c.m*j;
  c.level=(double *)R_alloc((size_t)c.width*3,sizeof(double));c.direction=c.level+c.width;c.adjoint=c.direction+c.width;
  c.output=PROTECT(allocVector(VECSXP,3));SET_VECTOR_ELT(c.output,0,allocMatrix(REALSXP,c.m,c.n));SET_VECTOR_ELT(c.output,1,allocMatrix(REALSXP,c.m,c.width));SET_VECTOR_ELT(c.output,2,ScalarReal(0));
  SEXP names=PROTECT(allocVector(STRSXP,3));SET_STRING_ELT(names,0,mkChar("weights"));SET_STRING_ELT(names,1,mkChar("projection"));SET_STRING_ELT(names,2,mkChar("ridge"));setAttrib(c.output,R_NamesSymbol,names);
  np_lp_full_row_workspace_init(&c.inverse);
  R_UnwindProtect(np_plot_alllarge_derivative_execute,&c,np_plot_alllarge_derivative_cleanup,&c,NULL);
  UNPROTECT(2);return c.output;
}
