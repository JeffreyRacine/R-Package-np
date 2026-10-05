#ifndef NP_CONDITIONAL_DELETED_QR_FAST_H
#define NP_CONDITIONAL_DELETED_QR_FAST_H
/* Production deleted-row QR helper. Caller-owned policy supplies lambda
 * and the pristine deleted Gram intercept; no rank rule lives here. */
typedef struct {
  int n,k,capacity,lwork;
  double *a,*scale,*tau,*v,*sqrtw,*work;
  int *pivot;
  /* Original donor identities for the conditional local-design factor. */
  int *row_index;
} NPConditionalQRFast;
static void np_cqr_fast_clear(NPConditionalQRFast *s) {
  free(s->a);free(s->scale);free(s->tau);free(s->v);free(s->sqrtw);
  free(s->work);free(s->pivot);free(s->row_index);memset(s,0,sizeof(*s));
}
static int np_cqr_fast_init(NPConditionalQRFast *s,int n,int k) {
  if(s->n==n && s->k==k && s->work)return 0;
  if(n<1||k<1||n>INT_MAX-k||k>(INT_MAX-1)/3)return 1;
  size_t m=(size_t)n+k;
  if(m>SIZE_MAX/(size_t)k/sizeof(double))return 1;
  np_cqr_fast_clear(s);s->n=n;s->k=k;s->capacity=(int)m;
  s->a=calloc(m*(size_t)k,sizeof(double));s->scale=calloc(k,sizeof(double));
  s->tau=calloc(k,sizeof(double));s->v=calloc(m,sizeof(double));
  s->sqrtw=calloc(n,sizeof(double));s->pivot=calloc(k,sizeof(int));
  s->row_index=calloc(n,sizeof(int));
  if(!s->a||!s->scale||!s->tau||!s->v||!s->sqrtw||!s->pivot||!s->row_index)goto fail;
  int lw=-1,info=0,one=1;double q1=0.,q2=0.,q3=0.;
  F77_CALL(dgeqp3)(&s->capacity,&k,s->a,&s->capacity,s->pivot,s->tau,&q1,&lw,&info);
  if(info)goto fail;
  F77_CALL(dormqr)("L","N",&s->capacity,&one,&k,s->a,&s->capacity,s->tau,
                   s->v,&s->capacity,&q2,&lw,&info FCONE FCONE);
  if(info)goto fail;
  F77_CALL(dgeqrf)(&s->capacity,&k,s->a,&s->capacity,s->tau,&q3,&lw,&info);
  if(info||!isfinite(q1)||!isfinite(q2)||!isfinite(q3)||q1>INT_MAX||q2>INT_MAX||q3>INT_MAX||q1<1||q2<1||q3<1)goto fail;
  s->lwork=(int)fmax(fmax(fmax(q1,q2),q3),3*k+1);
  s->work=calloc(s->lwork,sizeof(double));if(!s->work)goto fail;
  return 0;
fail:
  np_cqr_fast_clear(s);return 1;
}
static int np_cqr_fast_row(NPConditionalQRFast *s,int n,int k,const double *weights,
                          double **basis,int pos,double lambda,double anchor,
                          double *row) {
  if(pos<0||pos>=n||!isfinite(lambda)||lambda<0||
     (lambda>0&&(!isfinite(anchor)||anchor<=0))||np_cqr_fast_init(s,n,k))return 1;
  int m=lambda>0?s->capacity:n;
  if(m<k)return 1;
  double maxw=0.;
  for(int j=0;j<n;j++){
    if(!isfinite(weights[j])||weights[j]<0)return 1;
    if(j!=pos&&weights[j]>maxw)maxw=weights[j];
  }
  if(maxw==0)return 1;
  double mu=lambda/maxw;
  if(!isfinite(mu)||(lambda>0&&mu==0))return 1;
  double rootmu=sqrt(mu);
  for(int j=0;j<n;j++)s->sqrtw[j]=j==pos?0.:sqrt(weights[j]/maxw);
  int one=1,info=0;
  for(int l=0;l<k;l++){
    double *col=s->a+(size_t)s->capacity*l;
    for(int j=0;j<n;j++){
      if(!isfinite(basis[l][j]))return 1;
      col[j]=s->sqrtw[j]*basis[l][j];
    }
    if(lambda>0){memset(col+n,0,(size_t)k*sizeof(double));col[n+l]=rootmu;}
    s->scale[l]=F77_CALL(dnrm2)(&m,col,&one);
    if(!isfinite(s->scale[l])||s->scale[l]==0)return 1;
    for(int j=0;j<m;j++)col[j]/=s->scale[l];
    s->pivot[l]=0;
  }
  F77_CALL(dgeqp3)(&m,&k,s->a,&s->capacity,s->pivot,s->tau,s->work,&s->lwork,&info);
  if(info)return 1;
  memset(s->v,0,(size_t)m*sizeof(double));
  for(int i=0;i<k;i++){
    int p=s->pivot[i]-1;if(p<0||p>=k)return 1;
    double value=basis[p][pos]/s->scale[p];
    for(int j=0;j<i;j++)value-=s->a[j+(size_t)s->capacity*i]*s->v[j];
    double diagonal=s->a[i+(size_t)s->capacity*i];
    if(diagonal==0)return 1;
    s->v[i]=value/diagonal;if(!isfinite(s->v[i]))return 1;
  }
  F77_CALL(dormqr)("L","N",&m,&one,&k,s->a,&s->capacity,s->tau,s->v,&s->capacity,
                   s->work,&s->lwork,&info FCONE FCONE);
  if(info)return 1;
  double correction=0.;
  if(lambda>0){
    /* Last k rows encode sqrt(lambda/maxw) * u_normalized. */
    correction=(s->v[n]/rootmu)*(lambda/anchor);
    if(!isfinite(correction))return 1;
  }
  for(int j=0;j<n;j++){
    double value=s->sqrtw[j]*s->v[j];
    if(lambda>0&&j!=pos)value+=(weights[j]/maxw)*basis[0][j]*correction;
    if(!isfinite(value))return 1;
    row[j]=value;
  }
  return 0;
}

#endif
