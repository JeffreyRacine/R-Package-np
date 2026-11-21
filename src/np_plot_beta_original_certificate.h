/* Request-owned certificate for B' diag(w) B + rho I. The accepted solve,
 * its rank/ridge choice and the returned row are never changed here.
 * Include after np_plot_lp_precision_certificate.h (shared TwoSum primitive).
 * Expansions retain every term; exceptional products carry an outward radius.
 */
#include <fenv.h>
#include <stdint.h>

enum {
  NP_BETA_CERT_OK=1, NP_BETA_CERT_INACCURATE, NP_BETA_CERT_NONCONTRACTIVE,
  NP_BETA_CERT_STEP_LIMIT, NP_BETA_CERT_STATIONARY,
  NP_BETA_CERT_ARITHMETIC, NP_BETA_CERT_RESOURCE
};

typedef struct {
  double *e, radius;
  size_t length, capacity;
} NPBetaExpansion;
typedef struct {
  int n, p, error;
  size_t used, cells_count;
  double budget, q;
  NPBetaExpansion *cells, *gram, *z, *last, *residual, *forcing, *scratch;
  double *defect, *term, *radius, *next, *lower;
} NPBetaOriginalCertificate;

/* No persistent state: the small arithmetic check is made per request.
 * In particular, a build using unsafe reassociation cannot issue a certificate.
 */
static int np_beta_certificate_environment(void)
{
#if defined(__FAST_MATH__) || (defined(__FINITE_MATH_ONLY__) && __FINITE_MATH_ONLY__)
  return 0;
#else
  if(DBL_MANT_DIG!=53 || FLT_RADIX!=2 || FLT_EVAL_METHOD!=0 ||
     fegetround()!=FE_TONEAREST) return 0;
  volatile double one=1.0, small=0x1p-54, tiny=DBL_MIN;
  volatile double a=1.0+0x1p-52, b=1.0+0x1p-52, product=1.0+0x1p-51;
  double sum, residual;
  if(tiny*0.5==0.0 || !np_plot_two_sum(one,small,&sum,&residual) ||
     sum!=1.0 || residual!=0x1p-54 || fma(a,b,-product)!=0x1p-104 ||
     nextafter(one,INFINITY)-one!=0x1p-52 ||
     one-nextafter(one,0.0)!=0x1p-53) return 0;
  return 1;
#endif
}

static double np_beta_cert_finite(NPBetaOriginalCertificate *c,double x)
{
  if(!R_FINITE(x)) {c->error=NP_BETA_CERT_ARITHMETIC;return 0.0;}
  return x;
}
static double np_beta_cert_up(NPBetaOriginalCertificate *c,double x)
{return np_beta_cert_finite(c,nextafter(np_beta_cert_finite(c,x),INFINITY));}
static double np_beta_cert_add_up(NPBetaOriginalCertificate *c,double a,double b)
{return a==0.0?b:b==0.0?a:np_beta_cert_up(c,a+b);}
static double np_beta_cert_mul_up(NPBetaOriginalCertificate *c,double a,double b)
{return a==0.0||b==0.0?0.0:np_beta_cert_up(c,a*b);}
static double np_beta_cert_div_up(NPBetaOriginalCertificate *c,double a,double b)
{
  if(!(b>0.0)) {c->error=NP_BETA_CERT_ARITHMETIC;return 0.0;}
  return a==0.0?0.0:np_beta_cert_up(c,a/b);
}

static void *np_beta_cert_alloc(NPBetaOriginalCertificate *c,size_t n,size_t width)
{
  if(c->error)return NULL;
  if(width && n>SIZE_MAX/width) {c->error=NP_BETA_CERT_RESOURCE;return NULL;}
  size_t bytes=n*width;
  if(bytes>SIZE_MAX-c->used || (double)c->used+(double)bytes>c->budget) {
    c->error=NP_BETA_CERT_RESOURCE;return NULL;
  }
  void *p=calloc(n,width);
  if(p==NULL) {c->error=NP_BETA_CERT_RESOURCE;return NULL;}
  c->used+=bytes;return p;
}
static int np_beta_expansion_reserve(NPBetaOriginalCertificate *c,
                                    NPBetaExpansion *x,size_t count)
{
  if(c->error) return 0;
  if(count<=x->capacity) return 1;
  size_t cap=x->capacity?x->capacity:8;
  while(cap<count) {
    if(cap>SIZE_MAX/2) {c->error=NP_BETA_CERT_RESOURCE;return 0;}
    cap*=2;
  }
  if(cap>SIZE_MAX/sizeof(double)) {c->error=NP_BETA_CERT_RESOURCE;return 0;}
  const size_t extra=(cap-x->capacity)*sizeof(double);
  if(extra>SIZE_MAX-c->used || (double)c->used+(double)extra>c->budget) {
    c->error=NP_BETA_CERT_RESOURCE;return 0;
  }
  double *e=(double *)realloc(x->e,cap*sizeof(double));
  if(e==NULL) {c->error=NP_BETA_CERT_RESOURCE;return 0;}
  x->e=e;x->capacity=cap;c->used+=extra;return 1;
}
static void np_beta_expansion_value(NPBetaOriginalCertificate *c,
                                   NPBetaExpansion *x,double value)
{
  if(!np_beta_expansion_reserve(c,x,1)) return;
  x->e[0]=np_beta_cert_finite(c,value);x->length=1;x->radius=0.0;
}
static void np_beta_expansion_copy(NPBetaOriginalCertificate *c,
  NPBetaExpansion *out,const NPBetaExpansion *a)
{
  if(out==a || !np_beta_expansion_reserve(c,out,a->length)) return;
  memcpy(out->e,a->e,a->length*sizeof(double));
  out->length=a->length;out->radius=a->radius;
}
static void np_beta_expansion_grow(NPBetaOriginalCertificate *c,
                                  NPBetaExpansion *a,double value)
{
  if(a->length==SIZE_MAX || !np_beta_expansion_reserve(c,a,a->length+1)) return;
  double q=value;size_t count=0;
  for(size_t k=0;k<a->length;++k) {
    double sum,error;
    if(!np_plot_two_sum(q,a->e[k],&sum,&error)) {
      c->error=NP_BETA_CERT_ARITHMETIC;return;
    }
    if(error!=0.0) a->e[count++]=error;
    q=sum;
  }
  if(q!=0.0 || count==0) a->e[count++]=q;
  a->length=count;
}
/* out may alias a, never b. Each call site uses distinct scratch for products. */
static void np_beta_expansion_add(NPBetaOriginalCertificate *c,
  NPBetaExpansion *out,const NPBetaExpansion *a,const NPBetaExpansion *b,int sign)
{
  if(out==b) {c->error=NP_BETA_CERT_ARITHMETIC;return;}
  np_beta_expansion_copy(c,out,a);
  for(size_t k=0;k<b->length && !c->error;++k)
    np_beta_expansion_grow(c,out,sign*b->e[k]);
  out->radius=np_beta_cert_add_up(c,out->radius,b->radius);
}
static double np_beta_expansion_abs_sum(NPBetaOriginalCertificate *c,
                                       const NPBetaExpansion *a)
{
  double out=0.0;
  for(size_t k=0;k<a->length;++k) out=np_beta_cert_add_up(c,out,fabs(a->e[k]));
  return out;
}
static void np_beta_expansion_mul(NPBetaOriginalCertificate *c,
  NPBetaExpansion *out,const NPBetaExpansion *a,const NPBetaExpansion *b)
{
  if(out==a || out==b) {c->error=NP_BETA_CERT_ARITHMETIC;return;}
  np_beta_expansion_value(c,out,0.0);
  out->radius=np_beta_cert_add_up(c,np_beta_cert_add_up(c,
    np_beta_cert_mul_up(c,np_beta_expansion_abs_sum(c,a),b->radius),
    np_beta_cert_mul_up(c,np_beta_expansion_abs_sum(c,b),a->radius)),
    np_beta_cert_mul_up(c,a->radius,b->radius));
  for(size_t i=0;i<a->length && !c->error;++i)
    for(size_t j=0;j<b->length && !c->error;++j) {
      double x=a->e[i],y=b->e[j];if(x==0.0||y==0.0)continue;
      double product=np_beta_cert_finite(c,x*y);int ex,ey;
      (void)frexp(x,&ex);(void)frexp(y,&ey);
      if(ex+ey-2*DBL_MANT_DIG>=DBL_MIN_EXP-DBL_MANT_DIG) {
        double residual=np_beta_cert_finite(c,fma(x,y,-product));
        np_beta_expansion_grow(c,out,product);
        np_beta_expansion_grow(c,out,residual);
      } else {
        double lo=np_beta_cert_finite(c,nextafter(product,-INFINITY));
        double hi=np_beta_cert_finite(c,nextafter(product,INFINITY));
        out->radius=np_beta_cert_add_up(c,out->radius,fmax(product-lo,hi-product));
        np_beta_expansion_grow(c,out,product);
      }
    }
}
static void np_beta_expansion_scale(NPBetaOriginalCertificate *c,
  NPBetaExpansion *out,const NPBetaExpansion *a,double scale)
{
  NPBetaExpansion b={.e=&scale,.length=1,.radius=0.0};
  np_beta_expansion_mul(c,out,a,&b);
}
static NPPlotInterval np_beta_expansion_interval(NPBetaOriginalCertificate *c,
                                                const NPBetaExpansion *a)
{
  NPPlotInterval out={0.0,0.0};
  for(size_t k=0;k<a->length;++k) {
    double v=a->e[k];if(v==0.0)continue;
    if(out.lo==0.0&&out.hi==0.0) {out.lo=out.hi=v;continue;}
    out.lo=nextafter(np_beta_cert_finite(c,out.lo+v),-INFINITY);
    out.hi=nextafter(np_beta_cert_finite(c,out.hi+v),INFINITY);
  }
  if(a->radius!=0.0) {
    out.lo=nextafter(np_beta_cert_finite(c,out.lo-a->radius),-INFINITY);
    out.hi=nextafter(np_beta_cert_finite(c,out.hi+a->radius),INFINITY);
  }
  np_beta_cert_finite(c,out.lo);np_beta_cert_finite(c,out.hi);return out;
}
static double np_beta_expansion_abs(NPBetaOriginalCertificate *c,
                                    const NPBetaExpansion *a)
{return np_plot_iv_abs(np_beta_expansion_interval(c,a));}

static void np_beta_original_clear(NPBetaOriginalCertificate *c)
{
  if(c->cells)for(size_t i=0;i<c->cells_count;++i)free(c->cells[i].e);
  free(c->cells);free(c->defect);free(c->term);free(c->radius);
  free(c->next);free(c->lower);memset(c,0,sizeof(*c));
}
static int np_beta_original_init(NPBetaOriginalCertificate *c,int n,int p,double budget)
{
  c->n=n;c->p=p;c->budget=budget;
  if(!np_beta_certificate_environment()) return c->error=NP_BETA_CERT_ARITHMETIC;
  if(n<1||p<1||!R_FINITE(budget)||budget<=0 ||
     (size_t)p>(SIZE_MAX-8)/((size_t)p+4)) return c->error=NP_BETA_CERT_RESOURCE;
  c->cells_count=(size_t)p*p+4*(size_t)p+8;
  c->cells=np_beta_cert_alloc(c,c->cells_count,sizeof(*c->cells));
  c->defect=np_beta_cert_alloc(c,(size_t)p*p,sizeof(double));
  c->term=np_beta_cert_alloc(c,p,sizeof(double));
  c->radius=np_beta_cert_alloc(c,p,sizeof(double));
  c->next=np_beta_cert_alloc(c,p,sizeof(double));
  c->lower=np_beta_cert_alloc(c,n,sizeof(double));
  if(c->error)return c->error;
  c->gram=c->cells;c->z=c->gram+(size_t)p*p;c->last=c->z+p;
  c->residual=c->last+p;c->forcing=c->residual+p;c->scratch=c->forcing+p;
  return 0;
}

static void np_beta_expansion_dot(NPBetaOriginalCertificate *c,NPBetaExpansion *out,
  const NPBetaExpansion *a,const NPBetaExpansion *b)
{
  NPBetaExpansion *t=c->scratch;
  np_beta_expansion_value(c,out,0.0);
  for(int k=0;k<c->p && !c->error;++k) {
    np_beta_expansion_mul(c,t,a+k,b+k);np_beta_expansion_add(c,out,out,t,1);
  }
}
static void np_beta_expansion_ddot(NPBetaOriginalCertificate *c,NPBetaExpansion *out,
  const double *a,const NPBetaExpansion *b,size_t stride)
{
  NPBetaExpansion *t=c->scratch;
  np_beta_expansion_value(c,out,0.0);
  for(int k=0;k<c->p && !c->error;++k) {
    /* Keep the prototype's multiplication ordering: scalar first. */
    double v=a[k];NPBetaExpansion scalar={.e=&v,.length=1};
    np_beta_expansion_mul(c,t,&scalar,b+(size_t)k*stride);
    np_beta_expansion_add(c,out,out,t,1);
  }
}
static double np_beta_positive_dot(NPBetaOriginalCertificate *c,
                                  const double *a,const double *b)
{
  double s=0.0;
  for(int k=0;k<c->p;++k)
    s=np_beta_cert_add_up(c,s,np_beta_cert_mul_up(c,a[k],b[k]));
  return s;
}

/* H is the computed transposed inverse, stored by rows. It is a preconditioner,
 * never assumed exact. Retain Gram/E/q within this request for each direction.
 */
static int np_beta_original_prepare(NPBetaOriginalCertificate *c,
  double *const *B,const double *weights,double ridge,const double *H)
{
  if(c->error)return c->error;
  if(!np_beta_certificate_environment())return c->error=NP_BETA_CERT_ARITHMETIC;
  const int n=c->n,p=c->p;
  NPBetaExpansion *t=c->scratch,*v=t+1,*dot=t+2,*defect=t+3;
  for(int i=0;i<n;++i)
    if(B[0][i]!=1.0 || !R_FINITE(weights[i]) || weights[i]<0.0)
      return c->error=NP_BETA_CERT_ARITHMETIC;
  for(int i=0;i<p;++i)for(int j=i;j<p;++j) {
    NPBetaExpansion *a=c->gram+(size_t)i*p+j;
    np_beta_expansion_value(c,a,i==j?ridge:0.0);
    for(int k=0;k<n && !c->error;++k) {
      double w=weights[k];NPBetaExpansion scalar={.e=&w,.length=1};
      np_beta_expansion_scale(c,t,&scalar,B[i][k]);
      np_beta_expansion_scale(c,v,t,B[j][k]);np_beta_expansion_add(c,a,a,v,1);
    }
    if(i!=j)np_beta_expansion_copy(c,c->gram+(size_t)j*p+i,a);
  }
  c->q=0.0;
  for(int i=0;i<p && !c->error;++i) {
    double sum=0.0;
    for(int j=0;j<p;++j) {
      np_beta_expansion_ddot(c,dot,H+(size_t)i*p,c->gram+j,p);
      np_beta_expansion_value(c,defect,i==j?1.0:0.0);
      np_beta_expansion_add(c,defect,defect,dot,-1);
      double e=np_beta_expansion_abs(c,defect);c->defect[(size_t)i*p+j]=e;
      sum=np_beta_cert_add_up(c,sum,np_beta_cert_mul_up(c,e,1.0));
    }
    if(sum>c->q)c->q=sum;
  }
  return c->error?c->error:(c->q>=1.0?NP_BETA_CERT_NONCONTRACTIVE:0);
}

static void np_beta_original_map(NPBetaOriginalCertificate *c,
  double *const *B,const double *w,double ridge,double anchor,const double *returned,
  double *bound)
{
  NPBetaExpansion *t=c->scratch,*dot=t+1,*v=t+2,*scaled=t+3,*a=t+4,*b=t+5;
  double correction=np_beta_cert_add_up(c,1.0,np_beta_cert_div_up(c,ridge,anchor));
  for(int i=0;i<c->n && !c->error;++i) {
    if(w[i]==0.0) {bound[i]=c->lower[i]=fabs(returned[i]);continue;}
    np_beta_expansion_value(c,dot,0.0);
    for(int k=0;k<c->p;++k) {
      double value=B[k][i];NPBetaExpansion scalar={.e=&value,.length=1};
      np_beta_expansion_mul(c,t,&scalar,c->z+k);
      np_beta_expansion_add(c,dot,dot,t,1);
    }
    np_beta_expansion_scale(c,v,dot,w[i]);
    np_beta_expansion_value(c,a,returned[i]);np_beta_expansion_add(c,v,v,a,-1);
    if(ridge!=0.0) {
      np_beta_expansion_scale(c,scaled,v,anchor);
      np_beta_expansion_value(c,a,w[i]);np_beta_expansion_scale(c,b,a,B[0][i]);
      np_beta_expansion_scale(c,a,b,ridge);np_beta_expansion_mul(c,b,c->z,a);
      np_beta_expansion_add(c,scaled,scaled,b,1);v=scaled;
    }
    double delta=np_beta_expansion_abs(c,v);
    if(ridge!=0.0)delta=np_beta_cert_div_up(c,delta,anchor);
    NPPlotInterval interval=np_beta_expansion_interval(c,v);
    double low=interval.lo<=0.0&&interval.hi>=0.0?0.0:
      fmin(fabs(interval.lo),fabs(interval.hi));
    if(ridge!=0.0)low=fmax(0.0,nextafter(low/anchor,0.0));
    double propagated=0.0;
    for(int k=0;k<c->p;++k) {
      double r=np_beta_cert_mul_up(c,fabs(B[k][i]),
        np_beta_cert_mul_up(c,k==0?correction:1.0,c->radius[k]));
      propagated=np_beta_cert_add_up(c,propagated,np_beta_cert_mul_up(c,r,1.0));
    }
    propagated=np_beta_cert_mul_up(c,fabs(w[i]),propagated);
    bound[i]=np_beta_cert_add_up(c,delta,propagated);
    c->lower[i]=fmax(0.0,nextafter(low-propagated,-INFINITY));
    v=t+2;
  }
}

static double np_beta_original_ratio(NPBetaOriginalCertificate *c,
  const double *bound,const double *returned,int upper)
{
  double es=0.0,ms=0.0;
  for(int i=0;i<c->n;++i) {es=fmax(es,bound[i]);ms=fmax(ms,fabs(returned[i]));}
  if(es==0.0)return 0.0;
  if(ms==0.0)return upper?INFINITY:0.0;
  double se=0.0,sm=0.0;
  for(int i=0;i<c->n;++i) {
    if(upper) {
      double a=fmax(0.0,nextafter(fabs(returned[i])/ms,0.0));
      double b=bound[i]==0.0?0.0:np_beta_cert_up(c,bound[i]/es);
      sm=fmax(0.0,nextafter(sm+fmax(0.0,nextafter(a*a,0.0)),0.0));
      se=np_beta_cert_add_up(c,se,np_beta_cert_mul_up(c,b,b));
    } else {
      double a=fmax(0.0,nextafter(bound[i]/es,0.0));
      double b=returned[i]==0.0?0.0:np_beta_cert_up(c,fabs(returned[i])/ms);
      se=fmax(0.0,nextafter(se+fmax(0.0,nextafter(a*a,0.0)),0.0));
      sm=np_beta_cert_add_up(c,sm,np_beta_cert_mul_up(c,b,b));
    }
  }
  if(sm==0.0)return upper?INFINITY:0.0;
  if(upper)return np_beta_cert_mul_up(c,np_beta_cert_div_up(c,es,ms),
    np_beta_cert_up(c,sqrt(np_beta_cert_div_up(c,se,sm))));
  double ratio=fmax(0.0,nextafter(se/sm,0.0));
  double root=fmax(0.0,nextafter(sqrt(ratio),0.0));
  double scale=fmax(0.0,nextafter(es/ms,0.0));
  return fmax(0.0,nextafter(scale*root,0.0));
}

static int np_beta_original_evaluate(NPBetaOriginalCertificate *c,
  double *const *B,const double *w,double ridge,double anchor,const double *H,
  const double *direction,const double *initial,const double *returned,
  double *bounds,double *ratio,int *steps)
{
  *ratio=R_PosInf;*steps=0;
  if(c->error)return c->error;
  if(c->q>=1.0)return NP_BETA_CERT_NONCONTRACTIVE;
  const int p=c->p;
  double denominator=nextafter(1.0-c->q,-INFINITY);
  if(!(denominator>0.0)||!R_FINITE(anchor)||anchor<=0.0||!R_FINITE(ridge)||ridge<0.0)
    return c->error=NP_BETA_CERT_ARITHMETIC;
  for(int k=0;k<p;++k)np_beta_expansion_value(c,c->z+k,initial[k]);
  /* Initial private qualification ceiling; exhaustion is unavailable, never
   * a proof of inaccuracy or permission to change rank/ridge/the target. */
  for(int step=1;step<=32 && !c->error;++step) {
    NPBetaExpansion *dot=c->scratch+1,*value=c->scratch+2;
    for(int i=0;i<p;++i) {
      np_beta_expansion_dot(c,dot,c->gram+(size_t)i*p,c->z);
      np_beta_expansion_value(c,value,direction[i]);
      np_beta_expansion_add(c,c->residual+i,value,dot,-1);
    }
    for(int i=0;i<p;++i)
      np_beta_expansion_ddot(c,c->forcing+i,H+(size_t)i*p,c->residual,1);
    double forcing_norm=0.0,qpower=1.0;
    for(int i=0;i<p;++i) {
      np_beta_expansion_add(c,c->z+i,c->z+i,c->forcing+i,1);
      c->term[i]=np_beta_expansion_abs(c,c->forcing+i);c->radius[i]=0.0;
      forcing_norm=fmax(forcing_norm,c->term[i]);
    }
    for(int k=1;k<4;++k) {
      for(int i=0;i<p;++i)c->next[i]=np_beta_positive_dot(c,c->defect+(size_t)i*p,c->term);
      for(int i=0;i<p;++i) {
        c->term[i]=c->next[i];c->radius[i]=np_beta_cert_add_up(c,c->radius[i],c->term[i]);
      }
      qpower=np_beta_cert_mul_up(c,qpower,c->q);
    }
    double tail=np_beta_cert_div_up(c,np_beta_cert_mul_up(c,
      np_beta_cert_mul_up(c,qpower,c->q),forcing_norm),denominator);
    for(int i=0;i<p;++i)c->radius[i]=np_beta_cert_add_up(c,c->radius[i],tail);
    np_beta_original_map(c,B,w,ridge,anchor,returned,bounds);
    *ratio=np_beta_original_ratio(c,bounds,returned,1);
    double lower=np_beta_original_ratio(c,c->lower,returned,0);*steps=step;
    if(c->error)return c->error;
    if(lower>2e-11)return NP_BETA_CERT_INACCURATE;
    if(*ratio<=2e-11)return NP_BETA_CERT_OK;
    int same=step>1;
    for(int i=0;i<p&&same;++i) {
      NPBetaExpansion *a=c->last+i,*z=c->z+i;
      same=a->length==z->length&&a->radius==z->radius&&
        memcmp(a->e,z->e,z->length*sizeof(double))==0;
    }
    if(same)return NP_BETA_CERT_STATIONARY;
    for(int i=0;i<p;++i)np_beta_expansion_copy(c,c->last+i,c->z+i);
  }
  return c->error?c->error:NP_BETA_CERT_STEP_LIMIT;
}
