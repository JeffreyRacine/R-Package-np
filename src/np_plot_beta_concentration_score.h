/* Request-only inference arithmetic. No established beta owner uses this.
 * Stable concentration derivative of log Beta(u; 1+t*tau,1+(1-t)*tau).
 * Physical coordinates define the mathematical t/u; tau is the canonical
 * component concentration. Kernels and their signed coefficients stay owned
 * by the existing canonical beta primitives. */
#include <math.h>

static long double np_plot_beta_log1pmx(long double x)
{
  if(fabsl(x) > 0.125L) return log1pl(x)-x;
  long double term=x*x, out=-term/2.0L;
  for(int k=3;k<=32;++k) {
    term *= -x;
    const long double add=-term/(long double)k;
    const long double next=out+add;
    if(next==out) break;
    out=next;
  }
  return out;
}

/* psi(x)-log(x), x>=1. The first omitted asymptotic term after x^-20
 * is <1.5e-24 at x=16 and decreases thereafter. Recurrence is bounded
 * by 15 steps, avoiding the cancellation in psi(x)-log(x) at large x. */
static long double np_plot_beta_digamma_remainder(long double x)
{
  long double z=x, sum=0.0L;
  while(z<16.0L) {sum+=1.0L/z;z+=1.0L;}
  const long double v=1.0L/z,w=v*v;
  long double polynomial=174611.0L/6600.0L;
  polynomial=-43867.0L/14364.0L+w*polynomial;
  polynomial=3617.0L/8160.0L+w*polynomial;
  polynomial=-1.0L/12.0L+w*polynomial;
  polynomial=691.0L/32760.0L+w*polynomial;
  polynomial=-1.0L/132.0L+w*polynomial;
  polynomial=1.0L/240.0L+w*polynomial;
  polynomial=-1.0L/252.0L+w*polynomial;
  polynomial=1.0L/120.0L+w*polynomial;
  polynomial=-1.0L/12.0L+w*polynomial;
  return -v/2.0L+w*polynomial+log1pl((z-x)/x)-sum;
}

/* The units entry preserves a separately computed small beta tail when its
 * complement rounded to one (the CDF-score quadrature needs both tails). */
static double np_plot_beta_concentration_score_core(long double t,long double tc,
  long double u,long double uc,long double delta,long double tau)
{
  if(t<0.0L||tc<0.0L||u<=0.0L||uc<=0.0L||tau<0.0L) return NAN;
  if(t==0.0L) return (double)((u<=uc?log1pl(-u):logl(uc))+1.0L/(1.0L+tau));
  if(tc==0.0L) return (double)((uc<=u?log1pl(-uc):logl(u))+1.0L/(1.0L+tau));
  const long double dl=delta/t,dr=-delta/tc;
  long double entropy;
  if(fabsl(dl)<=0.125L&&fabsl(dr)<=0.125L)
    entropy=t*np_plot_beta_log1pmx(dl)+tc*np_plot_beta_log1pmx(dr);
  else {
    const long double lr1=fabsl(dl)<0.5L?log1pl(dl):logl(u)-logl(t);
    const long double lr2=fabsl(dr)<0.5L?log1pl(dr):logl(uc)-logl(tc);
    entropy=t*lr1+tc*lr2;
  }
  const long double a=1.0L+t*tau,b=1.0L+tc*tau,total=a+b;
  const long double ca=(2.0L*t-1.0L)/a,cb=(2.0L*tc-1.0L)/b;
  const long double la=fabsl(ca)<0.5L?log1pl(ca):logl(t)+logl(total)-logl(a);
  const long double lb=fabsl(cb)<0.5L?log1pl(cb):logl(tc)+logl(total)-logl(b);
  const long double location=t*la+tc*lb;
  const long double remainder=np_plot_beta_digamma_remainder(total)-
    t*np_plot_beta_digamma_remainder(a)-tc*np_plot_beta_digamma_remainder(b);
  return (double)(entropy+location+remainder);
}

static double np_plot_beta_concentration_score_units(double target,
  double target_complement,double observation,double observation_complement,
  double concentration)
{
  long double t=target,tc=target_complement,u=observation,uc=observation_complement;
  if(t<=tc)tc=1.0L-t;else t=1.0L-tc;
  if(u<=uc)uc=1.0L-u;else u=1.0L-uc;
  const long double delta=t<=tc?u-t:tc-uc;
  return np_plot_beta_concentration_score_core(t,tc,u,uc,delta,concentration);
}

static double np_plot_beta_concentration_score(double evaluation,
  double observation,double lower,double upper,double concentration)
{
  const long double q=evaluation,y=observation,L=(long double)upper-lower;
  const long double left=q-lower,right=(long double)upper-q;
  const long double yl=y-lower,yr=(long double)upper-y;
  if(!(L>0.0L)||left<0.0L||right<0.0L||yl<=0.0L||yr<=0.0L) return NAN;
  long double t,tc,u,uc;
  if(left<=right) {t=left/L;tc=1.0L-t;} else {tc=right/L;t=1.0L-tc;}
  if(yl<=yr) {u=yl/L;uc=1.0L-u;} else {uc=yr/L;u=1.0L-uc;}
  return np_plot_beta_concentration_score_core(t,tc,u,uc,(y-q)/L,concentration);
}
