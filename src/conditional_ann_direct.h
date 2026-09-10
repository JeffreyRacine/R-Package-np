/* Detached diagnostic: complete conditional influence, after fit totals.
 * A scaled Euclidean norm avoids squared-weight overflow and Gram subtraction.
 * No R/MPI dependency, allocation, clipping or alternate statistical target. */
#ifndef NP_CONDITIONAL_ANN_DIRECT_H
#define NP_CONDITIONAL_ANN_DIRECT_H
#include <math.h>
#include <stddef.h>

typedef struct {
  double scale, sumsq;
  int failed;
} NPANNConditionalNorm;

static inline int np_ann_direct_append(NPANNConditionalNorm *s,double x)
{
  double q;
  if(s->failed) return 0;
  x=fabs(x);
  if(!isfinite(x)) { s->failed=1; return 0; }
  if(x==0.0) return 1;
  if(s->scale<x) {
    q=s->scale/x;
    s->sumsq=1.0+s->sumsq*q*q;
    s->scale=x;
  } else {
    q=x/s->scale;
    s->sumsq+=q*q;
  }
  if(!isfinite(s->sumsq)) { s->failed=1;return 0; }
  return 1;
}

static inline double np_ann_direct_derivative(
  double a,double b,double residual,double t,double g)
{
  return (b-a*t)*residual-a*g;
}

static inline int np_ann_direct_finish(
  const NPANNConditionalNorm *s,size_t n,double *out)
{
  if(s->failed) return 0;
  if(n<=1) { *out=0.0; return 1; }
  *out=s->scale*sqrt(s->sumsq*((double)n/(double)(n-1)));
  return isfinite(*out);
}
#endif
