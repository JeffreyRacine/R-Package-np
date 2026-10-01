#ifndef NP_CONDITIONAL_DELETED_QR_H
#define NP_CONDITIONAL_DELETED_QR_H
/* Production deleted-row QR component with invocation-owned storage. */
/* Compensated reconstruction is selected by np_cqr_deleted_row; lambda and
 * the pristine deleted Gram intercept remain caller-owned policy. */
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include "conditional_deleted_qr_fast.h"
typedef struct { double hi,lo; } NPConditionalQRAcc;
static void np_cqr_add(NPConditionalQRAcc *s,double x){
  double t=s->hi+x,z=t-s->hi;
  double e=(s->hi-(t-z))+(x-z),v=s->lo+e,h=t+v;
  s->lo=v-(h-t);s->hi=h;
}
static void np_cqr_product(NPConditionalQRAcc *s,double x,double y){
  double v=x*y;np_cqr_add(s,v);np_cqr_add(s,fma(x,y,-v));
}

#include "conditional_deleted_dot.h"
typedef struct {
  NPConditionalQRFast fast;
  int m,k;
  double *q,*r,*scale,*v,*sw;
  int *pivot;
} NPConditionalQRDeleted;
static void np_cqr_deleted_mgs_clear(NPConditionalQRDeleted *s){
  free(s->q);free(s->r);free(s->scale);free(s->v);free(s->sw);free(s->pivot);
  s->q=s->r=s->scale=s->v=s->sw=NULL;s->pivot=NULL;s->m=s->k=0;
}
static void np_cqr_deleted_clear(NPConditionalQRDeleted *s){
  np_cqr_fast_clear(&s->fast);np_cqr_deleted_mgs_clear(s);
}
static int np_cqr_deleted_reserve(NPConditionalQRDeleted *s,int n,int k){
  if(n<1||k<1||n>INT_MAX-k)return 1;
  int m=n+k;
  if(s->m==m&&s->k==k&&s->q)return 0;
  if((size_t)m>SIZE_MAX/(size_t)k/sizeof(double))return 1;
  np_cqr_deleted_mgs_clear(s);
  s->q=calloc((size_t)m*k,sizeof(double));s->r=calloc((size_t)k*k,sizeof(double));
  s->scale=calloc(k,sizeof(double));s->v=calloc(k,sizeof(double));
  s->sw=calloc(n,sizeof(double));s->pivot=calloc(k,sizeof(int));
  if(!s->q||!s->r||!s->scale||!s->v||!s->sw||!s->pivot){np_cqr_deleted_mgs_clear(s);return 1;}
  s->m=m;s->k=k;return 0;
}
static int np_cqr_deleted_compensated(NPConditionalQRDeleted *s,int n,int k,
  const double *w,double **basis,int pos,double lambda,double anchor,double *row){
  if(np_cqr_deleted_reserve(s,n,k))return 1;
  const int m=lambda>0?s->m:n,one=1;double maxw=0;
  for(int i=0;i<n;i++)if(i!=pos&&w[i]>maxw)maxw=w[i];
  if(maxw<=0)return 1;
  double mu=lambda/maxw;
  if(!isfinite(mu)||(lambda>0&&mu==0))return 1;
  double *q=s->q,*r=s->r,*scale=s->scale,*v=s->v,*sw=s->sw;int *p=s->pivot;
  memset(r,0,(size_t)k*k*sizeof(double));
  for(int i=0;i<n;i++)sw[i]=i==pos?0:sqrt(w[i]/maxw);
  for(int a=0;a<k;a++){
    double *col=q+(size_t)m*a;
    for(int i=0;i<n;i++)col[i]=sw[i]*basis[a][i];
    if(lambda>0){memset(col+n,0,(size_t)k*sizeof(double));col[n+a]=sqrt(mu);}
    scale[a]=F77_CALL(dnrm2)(&m,col,&one);
    if(isfinite(scale[a])&&scale[a]>0)scale[a]=scalbn(1.,ilogb(scale[a]));
    if(!isfinite(scale[a])||scale[a]==0)return 1;
    for(int i=0;i<m;i++)col[i]/=scale[a];p[a]=a;
  }
  for(int a=0;a<k;a++){
    int pivot=a;double best=-1;
    for(int j=a;j<k;j++){
      double value=F77_CALL(dnrm2)(&m,q+(size_t)m*j,&one);
      if(!isfinite(value))return 1;
      if(value>best){best=value;pivot=j;}
    }
    if(pivot!=a){
      for(int i=0;i<m;i++){double t=q[i+(size_t)m*a];q[i+(size_t)m*a]=q[i+(size_t)m*pivot];q[i+(size_t)m*pivot]=t;}
      for(int j=0;j<a;j++){double t=r[j+(size_t)k*a];r[j+(size_t)k*a]=r[j+(size_t)k*pivot];r[j+(size_t)k*pivot]=t;}
      int t=p[a];p[a]=p[pivot];p[pivot]=t;
    }
    for(int j=0;j<a;j++){
      double z=np_cqr_dot4(q+(size_t)m*j,1,q+(size_t)m*a,1,m);r[j+(size_t)k*a]+=z;
      for(int i=0;i<m;i++)q[i+(size_t)m*a]=fma(-z,q[i+(size_t)m*j],q[i+(size_t)m*a]);
    }
    double norm=F77_CALL(dnrm2)(&m,q+(size_t)m*a,&one);
    if(!isfinite(norm)||norm==0)return 1;
    r[a+(size_t)k*a]=norm;for(int i=0;i<m;i++)q[i+(size_t)m*a]/=norm;
    for(int j=a+1;j<k;j++){
      double z=np_cqr_dot4(q+(size_t)m*a,1,q+(size_t)m*j,1,m);r[a+(size_t)k*j]=z;
      for(int i=0;i<m;i++)q[i+(size_t)m*j]=fma(-z,q[i+(size_t)m*a],q[i+(size_t)m*j]);
    }
  }
  for(int a=0;a<k;a++){
    NPConditionalQRAcc sum={basis[p[a]][pos]/scale[p[a]],0};
    for(int c=0;c<a;c++)np_cqr_product(&sum,-r[c+(size_t)k*a],v[c]);
    v[a]=(sum.hi+sum.lo)/r[a+(size_t)k*a];if(!isfinite(v[a]))return 1;
  }
  double correction=0;
  if(lambda>0)correction=(np_cqr_dot4(q+n,m,v,1,k)/sqrt(mu))*(lambda/anchor);
  if(!isfinite(correction))return 1;
  for(int i=0;i<n;i++){
    double z=sw[i]*np_cqr_dot4(q+i,m,v,1,k);
    if(lambda>0&&i!=pos)z+=(w[i]/maxw)*basis[0][i]*correction;
    if(!isfinite(z))return 1;
    row[i]=i==pos?0:z;
  }
  return 0;
}
static int np_cqr_deleted_row(NPConditionalQRDeleted *s,int n,int k,
  const double *w,double **basis,int pos,double lambda,double anchor,double *row,
  double *norm_out,int *owner_out){
  if(np_cqr_fast_row(&s->fast,n,k,w,basis,pos,lambda,anchor,row))return 1;
  int m=lambda>0?s->fast.capacity:n,one=1;
  double norm=F77_CALL(dnrm2)(&m,s->fast.v,&one);
  if(!isfinite(norm))return 1;
  if(norm_out)*norm_out=norm;if(owner_out)*owner_out=norm>1.;
  return norm>1.?np_cqr_deleted_compensated(s,n,k,w,basis,pos,lambda,anchor,row):0;
}

#endif
