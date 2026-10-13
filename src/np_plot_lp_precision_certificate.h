/* Scratch, request-only posterior numerical certificate. No fitted owner,
 * numerical policy, or accepted factor is modified. IEEE round-to-nearest
 * arithmetic without fast-math is required. nextafter gives outward bounds
 * for every elementary operation, including subnormal results. */
#include <math.h>
#include <float.h>

typedef struct { double lo,hi; } NPPlotInterval;
typedef struct {
  int p;
  double *inverse,*defect,*raw,*rlo,*rhi,*forcing,*error,*term,*next,*expansion;
} NPPlotLPCertificate;

static NPPlotInterval np_plot_iv(double x)
{NPPlotInterval a={x,x};return a;}
static NPPlotInterval np_plot_iv_add(NPPlotInterval a,NPPlotInterval b)
{
  if(a.lo==0.0&&a.hi==0.0)return b;
  if(b.lo==0.0&&b.hi==0.0)return a;
  NPPlotInterval c={nextafter(a.lo+b.lo,-INFINITY),nextafter(a.hi+b.hi,INFINITY)};
  return c;
}
static NPPlotInterval np_plot_iv_neg(NPPlotInterval a)
{NPPlotInterval c={-a.hi,-a.lo};return c;}
static NPPlotInterval np_plot_iv_mul(NPPlotInterval a,NPPlotInterval b)
{
  if((a.lo==0.0&&a.hi==0.0)||(b.lo==0.0&&b.hi==0.0))return np_plot_iv(0.0);
  double v[4]={a.lo*b.lo,a.lo*b.hi,a.hi*b.lo,a.hi*b.hi},lo=v[0],hi=v[0];
  for(int k=1;k<4;++k){if(v[k]<lo)lo=v[k];if(v[k]>hi)hi=v[k];}
  NPPlotInterval c={nextafter(lo,-INFINITY),nextafter(hi,INFINITY)};return c;
}
static double np_plot_iv_abs(NPPlotInterval a)
{return fmax(fabs(a.lo),fabs(a.hi));}
static double np_plot_up_add(double x,double y)
{if(x==0.0)return y;if(y==0.0)return x;return nextafter(x+y,INFINITY);}
static double np_plot_up_mul(double x,double y)
{if(x==0.0||y==0.0)return 0.0;return nextafter(x*y,INFINITY);}


/* Exact expansion of sums of binary64 products. TwoSum retains every bit
 * when all its intermediate values are finite (gradual underflow assumed).
 * The FMA residual is exact when the product's least possible binary quantum
 * is representable: frexp exponents ea+eb-2*DBL_MANT_DIG >= -1074.
 * We check this sufficient condition before TwoProd. Outside that domain,
 * fall back to the original per-operation outward interval; no capability or
 * numerical tolerance is silently substituted. The expansion workspace is
 * request-owned and reused across dots, never retained by an accepted fit. */
static int np_plot_two_sum(double a,double b,double *sum,double *error)
{
  const double s=a+b;if(!R_FINITE(s))return 0;
  const double bv=s-a, av=s-bv, br=b-bv, ar=a-av;
  const double e=ar+br;
  if(!R_FINITE(bv)||!R_FINITE(av)||!R_FINITE(br)||!R_FINITE(ar)||!R_FINITE(e))return 0;
  *sum=s;*error=e;return 1;
}
static int np_plot_expansion_grow(double *e,int *length,double value)
{
  double q=value;int out=0;
  for(int k=0;k<*length;++k){double s,r;
    if(!np_plot_two_sum(q,e[k],&s,&r))return 0;
    if(r!=0.0)e[out++]=r;
    q=s;
  }
  if(q!=0.0||out==0)e[out++]=q;
  *length=out;return 1;
}
static int np_plot_product_expansion(double *e,int *length,double a,double b)
{
  if(!R_FINITE(a)||!R_FINITE(b))return 0;
  if(a==0.0||b==0.0)return 1;
  int ea,eb;(void)frexp(a,&ea);(void)frexp(b,&eb);
  /* Require the entire exact product, including its residual, to have a
   * representable least bit. No reliance on platform long-double width. */
  if(ea+eb-2*DBL_MANT_DIG < DBL_MIN_EXP-DBL_MANT_DIG)return 0;
  const double product=a*b;if(!R_FINITE(product))return 0;
  const double residual=fma(a,b,-product);if(!R_FINITE(residual))return 0;
  return np_plot_expansion_grow(e,length,product)&&
         np_plot_expansion_grow(e,length,residual);
}
/* Return an enclosure of initial-sum(a[k*sa]*b[k*sb]). Successful EFT
 * produces an exact expansion; only its final, normally short sum is rounded
 * outward. Both the tight and fallback paths target the exact same dot. */
static NPPlotInterval np_plot_lp_exact_residual_dot(double *expansion,int n,
  double initial,const double *a,size_t sa,const double *b,size_t sb)
{
  int length=0,ok=R_FINITE(initial);
  if(ok)ok=np_plot_expansion_grow(expansion,&length,initial);
  for(int k=0;k<n&&ok;++k)
    ok=np_plot_product_expansion(expansion,&length,-a[(size_t)k*sa],b[(size_t)k*sb]);
  NPPlotInterval result=np_plot_iv(0.0);
  if(ok){
    for(int k=0;k<length;++k)result=np_plot_iv_add(result,np_plot_iv(expansion[k]));
    return result;
  }
  result=np_plot_iv(initial);
  for(int k=0;k<n;++k)result=np_plot_iv_add(result,np_plot_iv_neg(np_plot_iv_mul(
    np_plot_iv(a[(size_t)k*sa]),np_plot_iv(b[(size_t)k*sb]))));
  return result;
}

/* Allocate once per request, not once per evaluation row. The enclosing
 * request budget must include 16*p*p+96*p+32 bytes for this workspace. */
static void np_plot_lp_certificate_init(NPPlotLPCertificate *c,int p)
{
  c->p=p;c->inverse=(double *)R_alloc((size_t)p*p,sizeof(double));
  c->defect=(double *)R_alloc((size_t)p*p,sizeof(double));
  c->raw=(double *)R_alloc(p,sizeof(double));c->rlo=(double *)R_alloc(p,sizeof(double));
  c->rhi=(double *)R_alloc(p,sizeof(double));c->forcing=(double *)R_alloc(p,sizeof(double));
  c->error=(double *)R_alloc(p,sizeof(double));c->term=(double *)R_alloc(p,sizeof(double));
  c->next=(double *)R_alloc(p,sizeof(double));
  c->expansion=(double *)R_alloc((size_t)2*p+2,sizeof(double));
}

/* Internal arithmetic stage, also callable by a standalone independent
 * certificate oracle. G is the exact retained POST-ridge Gram in column
 * major storage; inverse is its computed inverse from the retained factor.
 * raw is an uncorrected computed G^{-T} d. Neither is assumed exact. */
static int np_plot_lp_certificate_arithmetic(NPPlotLPCertificate *c,
  const double *G,const double *d,double *defect_norm)
{
  const int p=c->p;double q=0.0;
  /* E=I-H*A, H=inverse^T, A=G^T. Store |E| upper bounds. */
  for(int i=0;i<p;++i){
    double sum=0.0;
    for(int j=0;j<p;++j){
      NPPlotInterval v=np_plot_lp_exact_residual_dot(c->expansion,p,
        i==j?1.0:0.0,c->inverse+(size_t)i*p,1,G+j,(size_t)p);
      const double t=np_plot_iv_abs(v);if(!R_FINITE(t))return 0;
      c->defect[i+(size_t)j*p]=t;sum=np_plot_up_add(sum,t);
    }
    if(sum>q)q=sum;
  }
  *defect_norm=q;if(!R_FINITE(q)||q>=1.0)return 0;
  /* Residual enclosure r=d-G^T raw includes arithmetic uncertainty. */
  for(int i=0;i<p;++i){
    NPPlotInterval r=np_plot_lp_exact_residual_dot(c->expansion,p,d[i],
      G+(size_t)i*p,1,c->raw,1);
    if(!R_FINITE(r.lo)||!R_FINITE(r.hi))return 0;
    c->rlo[i]=r.lo;c->rhi[i]=r.hi;
  }
  double forcing_norm=0.0;
  for(int i=0;i<p;++i){
    NPPlotInterval f=np_plot_iv(0.0);
    for(int j=0;j<p;++j){NPPlotInterval r={c->rlo[j],c->rhi[j]};
      f=np_plot_iv_add(f,np_plot_iv_mul(np_plot_iv(c->inverse[j+(size_t)i*p]),r));}
    const double a=np_plot_iv_abs(f);if(!R_FINITE(a))return 0;
    c->forcing[i]=c->error[i]=c->term[i]=a;if(a>forcing_norm)forcing_norm=a;
  }
  /* Four nonnegative Neumann terms plus a proved infinity-norm tail. */
  double qpower=1.0;
  for(int step=1;step<4;++step){
    for(int i=0;i<p;++i){double s=0.0;
      for(int j=0;j<p;++j)s=np_plot_up_add(s,np_plot_up_mul(c->defect[i+(size_t)j*p],c->term[j]));
      c->next[i]=s;}
    for(int i=0;i<p;++i){c->term[i]=c->next[i];c->error[i]=np_plot_up_add(c->error[i],c->next[i]);}
    qpower=np_plot_up_mul(qpower,q);
  }
  qpower=np_plot_up_mul(qpower,q);
  const double denominator=nextafter(1.0-q,-INFINITY);
  if(!(denominator>0.0))return 0;
  const double tail=nextafter(np_plot_up_mul(qpower,forcing_norm)/denominator,INFINITY);
  if(!R_FINITE(tail))return 0;
  for(int i=0;i<p;++i){c->error[i]=np_plot_up_add(c->error[i],tail);if(!R_FINITE(c->error[i]))return 0;}
  return 1;
}

/* Certify one accepted level/derivative adjoint without refactorization.
 * accepted is the already returned canonical CORRECTED adjoint; bound gets
 * absolute component errors relative to the exact frozen-geometry operator.
 * This mutates only RHS scratch. gram_source and the retained LU are kept. */
static int np_plot_lp_certificate_adjoint(NPPlotLPCertificate *c,
  NPLPSolveWorkspace *w,const double *d,const double *accepted,double anchor,
  const NPLPSolvePolicyDiagnostics *policy,double *bound,double *defect_norm)
{
  const int p=c->p;
  for(int i=0;i<p;++i)w->rhs_source[i]=d[i];
  if(!np_lp_solve_workspace_solve_transpose_factored(w,p,1))return 0;
  for(int i=0;i<p;++i)c->raw[i]=w->rhs_work[i];
  for(int k=0;k<p;++k){
    for(int i=0;i<p;++i)w->rhs_source[i]=i==k?1.0:0.0;
    if(!np_lp_solve_workspace_solve_factored(w,p,1))return 0;
    for(int i=0;i<p;++i){double v=w->rhs_work[i];if(!R_FINITE(v))return 0;c->inverse[i+(size_t)k*p]=v;}
  }
  if(!np_plot_lp_certificate_arithmetic(c,w->gram_source,d,defect_norm))return 0;
  NPPlotInterval correction=np_plot_iv(1.0);
  if(policy->ridge_total!=0.0){
    if(!R_FINITE(anchor)||anchor==0.0)return 0;
    const double z=policy->ridge_total/anchor;
    NPPlotInterval ratio={nextafter(z,-INFINITY),nextafter(z,INFINITY)};
    correction=np_plot_iv_add(correction,ratio);
  }
  for(int i=0;i<p;++i){
    NPPlotInterval truth={nextafter(c->raw[i]-c->error[i],-INFINITY),nextafter(c->raw[i]+c->error[i],INFINITY)};
    if(i==0)truth=np_plot_iv_mul(correction,truth);
    truth=np_plot_iv_add(truth,np_plot_iv(-accepted[i]));
    bound[i]=np_plot_iv_abs(truth);if(!R_FINITE(bound[i]))return 0;
  }
  return 1;
}

/* Enclose K_i * B_i^T a with each a_k in accepted[k]+/-bound[k].
 * computed is the actual returned row entry, not a freshly rounded center.
 * Calling code can retain these elementwise bounds to propagate through
 * LSQ residual and full covariance factors; no unit absolute floor is used. */
static double np_plot_lp_certificate_map_error(int p,double **basis,int i,
  double weight,const double *accepted,const double *bound,double computed)
{
  NPPlotInterval v=np_plot_iv(0.0);
  for(int k=0;k<p;++k){NPPlotInterval a={nextafter(accepted[k]-bound[k],-INFINITY),nextafter(accepted[k]+bound[k],INFINITY)};
    v=np_plot_iv_add(v,np_plot_iv_mul(np_plot_iv(basis[k][i]),a));}
  v=np_plot_iv_mul(np_plot_iv(weight),v);
  return np_plot_iv_abs(np_plot_iv_add(v,np_plot_iv(-computed)));
}
