/* Signed GNN polynomial arithmetic. Incumbent rank/ridge admission is owned
 * by the caller. This equivalent zero-ridge representation never regularizes. */
typedef struct {
  int n, k;
  double *basis, *gram, *solution;
  int *pivot;
  NPConditionalQRAcc *polynomial;
} NPGNNSignedRow;
static void np_cgnn_signed_clear(NPGNNSignedRow *s) {
  free(s->basis); free(s->gram); free(s->solution); free(s->pivot);
  free(s->polynomial); memset(s,0,sizeof(*s));
}
static int np_cgnn_signed_prepare(NPGNNSignedRow *s,int n,int k) {
  size_t nk,kk,pk;
  if(n<1 || k<2 || !np_size_mul_checked(n,k,&nk) ||
     !np_size_mul_checked(k,k,&kk) || !np_size_mul_checked(k,3,&pk))return 1;
  s->n=n;s->k=k;
  s->basis=np_cgnn_calloc(nk,sizeof(double));
  s->gram=np_cgnn_calloc(kk,sizeof(double));
  s->solution=np_cgnn_calloc(k,sizeof(double));
  s->pivot=np_cgnn_calloc(k,sizeof(int));
  s->polynomial=np_cgnn_calloc(pk,sizeof(NPConditionalQRAcc));
  return !(s->basis && s->gram && s->solution && s->pivot && s->polynomial);
}
static int np_cgnn_signed_row(NPGNNSignedRow *s,const double *x,const double *w,
    int held,const int *terms,NPConditionalQRAcc *coef,double *row,double *radius) {
  const int n=s->n,k=s->k;double lo=x[held],hi=lo;
  for(int p=0;p<k;p++)if(terms[p]!=p)return 1;
  for(int j=0;j<n;j++)if(w[j]!=0){lo=fmin(lo,x[j]);hi=fmax(hi,x[j]);}
  const double center=.5*lo+.5*hi;
  const double scale=fmax(fabs(lo-center),fabs(hi-center));
  if(!(scale>0) || !R_FINITE(scale))return 1;
  for(int j=0;j<n;j++){
    const double u=(x[j]-center)/scale;
    s->basis[j]=1.;s->basis[j+(size_t)n]=u;
    for(int p=2;p<k;p++)s->basis[j+(size_t)n*p]=
      2*u*s->basis[j+(size_t)n*(p-1)]-s->basis[j+(size_t)n*(p-2)];
  }
  for(int p=0;p<k;p++){
    s->solution[p]=s->basis[held+(size_t)n*p];
    for(int q=0;q<=p;q++){
      NPConditionalQRAcc sum={0,0};
      for(int j=0;j<n;j++)if(w[j]!=0){
        NPConditionalQRAcc term=cqr_dd_mul((NPConditionalQRAcc){w[j],0},
          cqr_dd_mul((NPConditionalQRAcc){s->basis[j+(size_t)n*p],0},
                     (NPConditionalQRAcc){s->basis[j+(size_t)n*q],0}));
        np_cqr_add(&sum,term.hi);np_cqr_add(&sum,term.lo);
      }
      s->gram[p+(size_t)k*q]=s->gram[q+(size_t)k*p]=sum.hi+sum.lo;
    }
  }
  int nrhs=1,info=0;
  /* No second rank test or ridge retry: the incumbent decision is retained. */
  F77_CALL(dgesv)(&k,&nrhs,s->gram,&k,s->pivot,s->solution,&k,&info);
  if(info!=0)return 1;
  double denominator;
  NPConditionalQRAcc self={0,0};
  for(int p=0;p<k;p++)np_cqr_product(&self,s->basis[held+(size_t)n*p],s->solution[p]);
  if(!np_lp_delete_denominator((self.hi+self.lo)*w[held],&denominator))return 1;
  for(int j=0;j<n;j++){
    NPConditionalQRAcc sum={0,0};
    if(j!=held && w[j]!=0)for(int p=0;p<k;p++)
      np_cqr_product(&sum,s->basis[j+(size_t)n*p],s->solution[p]);
    row[j]=(sum.hi+sum.lo)*w[j]/denominator;
    if(!R_FINITE(row[j]))return 1;
  }
  /* Convert T_p(offset+u) to target-centered powers in twofold arithmetic.
   * The integral consumes exactly this signed influence polynomial. */
  NPConditionalQRAcc offset={x[held],0};np_cqr_add(&offset,-center);
  offset=cqr_dd_div(offset,scale);
  NPConditionalQRAcc *prev=s->polynomial,*cur=prev+k,*next=cur+k;
  memset(s->polynomial,0,(size_t)3*k*sizeof(*prev));
  memset(coef,0,(size_t)k*sizeof(*coef));prev[0]=(NPConditionalQRAcc){1,0};
  cur[0]=offset;cur[1]=(NPConditionalQRAcc){1,0};
  for(int p=0;p<k;p++){
    NPConditionalQRAcc *poly=p==0?prev:cur;
    const NPConditionalQRAcc weight=cqr_dd_div((NPConditionalQRAcc){s->solution[p],0},denominator);
    for(int q=0;q<=p;q++){
      NPConditionalQRAcc v=cqr_dd_mul(poly[q],weight);
      np_cqr_add(coef+q,v.hi);np_cqr_add(coef+q,v.lo);
    }
    if(p>0 && p+1<k){
      memset(next,0,(size_t)k*sizeof(*next));
      for(int q=0;q<=p;q++){
        NPConditionalQRAcc v=cqr_dd_mul(offset,cur[q]);
        np_cqr_add(next+q,2*v.hi);np_cqr_add(next+q,2*v.lo);
        np_cqr_add(next+q+1,2*cur[q].hi);np_cqr_add(next+q+1,2*cur[q].lo);
        np_cqr_add(next+q,-prev[q].hi);np_cqr_add(next+q,-prev[q].lo);
      }
      NPConditionalQRAcc *tmp=prev;prev=cur;cur=next;next=tmp;
    }
  }
  for(int p=0;p<k;p++)if(!R_FINITE(coef[p].hi)||!R_FINITE(coef[p].lo))return 1;
  *radius=scale;return 0;
}
