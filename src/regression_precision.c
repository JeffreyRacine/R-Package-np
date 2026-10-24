#include "regression_precision.h"
#include "regression_contrast.h"
#include <R.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

void np_regression_precision_init(NPRegressionPrecision *q)
{ memset(q, 0, sizeof(*q)); }
void np_regression_precision_clear(NPRegressionPrecision *q)
{
  free(q->a); free(q->vectors); free(q->work); free(q->pivot);
  np_regression_precision_init(q);
}

/* Conservative inverse-norm bounds for the retained triangular factors.
 * Permutations preserve the infinity norm. Scratch is already owned by the
 * accepted solve. This screen neither changes that factor nor adjudicates rank.
 * It is deliberately independent of the response and requested outputs. */
int np_regression_precision_needed(NPLPSolveWorkspace *s, int p)
{
  if(p <= 1) return 0;
  if(!s->factor_ready || s->factor_p != p || s->rank_values == NULL) return 1;
  double norm=0.0, lower=0.0, upper=0.0;
  double *v=s->rank_values;
  for(int i=0;i<p;++i) {
    double z=0.0;
    for(int j=0;j<p;++j) z+=fabs(s->gram_source[i+(size_t)p*j]);
    if(z>norm) norm=z;
    z=1.0;
    for(int j=0;j<i;++j) z+=fabs(s->gram_work[i+(size_t)p*j])*v[j];
    v[i]=z;if(z>lower)lower=z;
  }
  for(int i=p-1;i>=0;--i) {
    double z=1.0;
    for(int j=i+1;j<p;++j) z+=fabs(s->gram_work[i+(size_t)p*j])*v[j];
    z/=fabs(s->gram_work[i+(size_t)p*i]);v[i]=z;if(z>upper)upper=z;
  }
  const double gamma=(p*DBL_EPSILON)/(1.0-p*DBL_EPSILON);
  const double bound=gamma*norm*lower*upper;
  return !R_FINITE(bound) || bound>1e-10;
}

static int reserve_double(double **x, size_t *capacity, size_t count)
{
  if(count<=*capacity)return 1;
  if(count>SIZE_MAX/sizeof(double))return 0;
  double *next=(double *)realloc(*x,count*sizeof(double));
  if(next==NULL)return 0;
  *x=next;*capacity=count;return 1;
}
static int reserve_work(NPRegressionPrecision *q, double query)
{
  if(!R_FINITE(query)||query<1||query>INT_MAX)return 0;
  int need=(int)query;
  if(need>q->work_capacity) {
    double *next=(double *)realloc(q->work,(size_t)need*sizeof(double));
    if(next==NULL)return 0;
    q->work=next;q->work_capacity=need;
  }
  return 1;
}
int np_regression_precision_factor(NPRegressionPrecision *q,
  double *const *basis, const double *weights, int n, int p,
  double ridge, double anchor)
{
  q->weights=NULL;
  if(n<1||p<1||n>INT_MAX-p||!R_FINITE(ridge)||ridge<0||
     !R_FINITE(anchor)||anchor<=0||basis==NULL||weights==NULL)return 0;
  int m=ridge>0?n+p:n;
  if(m<p||(size_t)m>SIZE_MAX/(size_t)p)return 0;
  if(!reserve_double(&q->a,&q->matrix_capacity,(size_t)m*p)||
     !reserve_double(&q->vectors,&q->vector_capacity,2*(size_t)p+2*(size_t)n+m))return 0;
  if(p>q->pivot_capacity) {
    int *next=(int *)realloc(q->pivot,(size_t)p*sizeof(int));
    if(next==NULL)return 0;
    q->pivot=next;q->pivot_capacity=p;
  }
  q->n=n;q->p=p;q->m=m;
  q->tau=q->vectors;q->scale=q->tau+p;q->sw=q->scale+p;
  q->u=q->sw+n;q->row=q->u+m;
  double maxw=0.0;
  for(int i=0;i<n;++i){
    if(!R_FINITE(weights[i])||weights[i]<0)return 0;
    if(weights[i]>maxw)maxw=weights[i];
  }
  if(maxw==0)return 0;
  q->max_weight=maxw;q->anchor=anchor;q->root_ridge=sqrt(ridge/maxw);
  if(!R_FINITE(q->root_ridge))return 0;
  for(int i=0;i<n;++i)q->sw[i]=sqrt(weights[i]/maxw);
  for(int k=0;k<p;++k){
    double norm=0;
    for(int i=0;i<n;++i){
      if(!R_FINITE(basis[k][i]))return 0;
      double z=q->sw[i]*basis[k][i];q->a[i+(size_t)m*k]=z;
    }
    if(ridge>0){
      memset(q->a+n+(size_t)m*k,0,(size_t)p*sizeof(double));
      q->a[n+k+(size_t)m*k]=q->root_ridge;
    }
    const int one=1;norm=F77_CALL(dnrm2)(&m,q->a+(size_t)m*k,&one);
    if(norm==0||!R_FINITE(norm))return 0;
    q->scale[k]=norm;q->pivot[k]=0;
    for(int i=0;i<m;++i)q->a[i+(size_t)m*k]/=norm;
  }
  int info=0,lwork=-1;double query;
  const int query_workspace=q->work_m!=m || q->work_p!=p;
  if(query_workspace) {
    F77_CALL(dgeqp3)(&m,&p,q->a,&m,q->pivot,q->tau,&query,&lwork,&info);
    if(info||!reserve_work(q,query))return 0;
  }
  F77_CALL(dgeqp3)(&m,&p,q->a,&m,q->pivot,q->tau,q->work,&q->work_capacity,&info);
  if(info)return 0;
  int one=1;lwork=-1;
  if(query_workspace) {
    F77_CALL(dormqr)("L","N",&m,&one,&p,q->a,&m,q->tau,q->u,&m,
                    &query,&lwork,&info FCONE FCONE);
    if(info||!reserve_work(q,query))return 0;
    q->work_m=m;q->work_p=p;
  }
  q->weights=weights;
  return 1;
}
int np_regression_precision_row(NPRegressionPrecision *q, const double *direction)
{
  if(q->weights==NULL||direction==NULL)return 0;
  const int m=q->m,p=q->p,n=q->n;
  memset(q->u,0,(size_t)m*sizeof(double));
  for(int i=0;i<p;++i){
    int k=q->pivot[i]-1;
    double z=direction[k]/q->scale[k];
    for(int j=0;j<i;++j)z-=q->a[j+(size_t)m*i]*q->u[j];
    double diag=q->a[i+(size_t)m*i];
    if(diag==0||!R_FINITE(z))return 0;
    q->u[i]=z/diag;if(!R_FINITE(q->u[i]))return 0;
  }
  int one=1,info;
  F77_CALL(dormqr)("L","N",&m,&one,&p,q->a,&m,q->tau,q->u,&m,
                  q->work,&q->work_capacity,&info FCONE FCONE);
  if(info)return 0;
  const double correction=q->root_ridge>0?q->root_ridge*q->u[n]/q->anchor:0;
  for(int i=0;i<n;++i){
    q->row[i]=q->sw[i]*q->u[i]+correction*q->weights[i];
    if(!R_FINITE(q->row[i]))return 0;
  }
  return 1;
}
double np_regression_precision_apply(const NPRegressionPrecision *q,
  const double *response)
{
  NPContrastNumber sum={0};
  for(int i=0;i<q->n;++i)
    sum=np_contrast_add(sum,np_contrast_two_product(q->row[i],response[i]));
  return np_contrast_value(sum);
}
