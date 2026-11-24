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
  free(q->local_columns); free(q->a); free(q->vectors); free(q->work); free(q->pivot);
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
/* The four pure-cubic GLP columns depend on one predictor.
 * Three distinct supported column tuples prove one exact dependency. The
 * triangular checks only guard corroboration; they never select a fit rank. */
static int cubic_support_alias(NPRegressionPrecision *q, double *const *basis)
{
  if(q->p != 7 || q->root_ridge != 0.0) return 0;
  const int stride = q->m + 1;
  const double guard = q->m * DBL_EPSILON * fabs(q->a[0]);
  if(fabs(q->a[6 * stride]) > guard || fabs(q->a[5 * stride]) <= guard)
    return 0;

  /* The GLP lower set is enumerated by increasing second-coordinate degree. */
  for(int k = 0; k < 4; ++k)
    q->cubic_columns[k] = q->check_cubic_support == 1 ? k : 2*k;
  int count = 0;
  for(int i = 0; i < q->n; ++i) {
    if(q->weights[i] == 0.0) continue;
    int level = 0;
    while(level < count && q->cubic_levels[level][1] != basis[q->cubic_columns[1]][i]) ++level;
    if(level == count) {
      if(count == 3) return 0;
      for(int k = 0; k < 4; ++k) q->cubic_levels[count][k] = basis[q->cubic_columns[k]][i];
      ++count;
    } else {
      for(int k = 0; k < 4; ++k)
        if(q->cubic_levels[level][k] != basis[q->cubic_columns[k]][i]) return 0;
    }
  }
  return count == 3;
}

/* A target with no pure-cubic component, or evaluation at a supported
 * abscissa, annihilates the alias. A different functional retains the
 * incumbent public fit convention. Corrected bands check it separately. */
static int cubic_identified_direction(const NPRegressionPrecision *q,
                                      const double *v)
{
  const int *col = q->cubic_columns;
  if(v[col[0]] == 0.0 && v[col[1]] == 0.0 &&
     v[col[2]] == 0.0 && v[col[3]] == 0.0) return 1;
  for(int level = 0; level < 3; ++level) {
    int k = 0;
    while(k < 4 && q->cubic_levels[level][k] == v[col[k]]) ++k;
    if(k == 4) return 1;
  }
  return 0;
}

/* The first normalized Legendre coordinate is affine in physical x. Center
 * at the highest-weight donor and scale by the weighted RMS distance. This
 * changes coordinates only after the native full-rank zero-ridge decision. */
static int local_cubic_columns(NPRegressionPrecision *q, const double *linear,
                               const double *weights, int n, double maxw,
                               double **columns)
{
  if(!reserve_double(&q->local_columns,&q->local_capacity,(size_t)n*4))return 0;
  int peak=0;
  while(peak<n && weights[peak]!=maxw)++peak;
  if(peak==n)return 0;
  const double center=linear[peak];
  double mass=0, square=0;
  for(int i=0;i<n;++i) {
    const double w=weights[i]/maxw, dx=linear[i]-center;
    mass+=w; square+=w*dx*dx;
  }
  const double scale=sqrt(square/mass);
  if(!(scale>0) || !R_FINITE(scale))return 0;
  for(int k=0;k<4;++k)columns[k]=q->local_columns+(size_t)n*k;
  for(int i=0;i<n;++i) {
    const double t=weights[i]>0?(linear[i]-center)/scale:0;
    columns[0][i]=1; columns[1][i]=t;
    columns[2][i]=t*t; columns[3][i]=t*t*t;
  }
  q->cubic_center=center; q->cubic_scale=scale;
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
  double *local[4];
  q->cubic_scale=0;
  if(q->local_cubic && ridge==0.0) {
    if(p!=4 || !local_cubic_columns(q,basis[1],weights,n,maxw,local))return 0;
    basis=local;
  }
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
  q->cubic_alias = q->check_cubic_support && cubic_support_alias(q,basis);
  return 1;
}
int np_regression_precision_row(NPRegressionPrecision *q, const double *direction)
{
  if(q->weights==NULL||direction==NULL)return 0;
  double local[4];
  if(q->cubic_scale>0) {
    /* In the normalized Legendre basis B2=sqrt(5)*(B1^2-1)/2.
     * A first-derivative request thus supplies B1(x) through B2'/B1'. */
    const double scale=q->cubic_scale;
    if(direction[0]==1.0) {
      const double t=(direction[1]-q->cubic_center)/scale;
      local[0]=1; local[1]=t; local[2]=t*t; local[3]=t*t*t;
    } else if(direction[0]==0.0 && direction[1]!=0.0) {
      const double t=(direction[2]/(sqrt(5.0)*direction[1])-q->cubic_center)/scale;
      const double slope=direction[1]/scale;
      local[0]=0; local[1]=slope; local[2]=2*t*slope; local[3]=3*t*t*slope;
    } else return 0;
    direction=local;
  }
  const int m=q->m,n=q->n;
  const int p = q->cubic_alias && cubic_identified_direction(q,direction)
    ? q->p - 1 : q->p;
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
