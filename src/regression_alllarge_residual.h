#ifndef NP_REGRESSION_ALLLARGE_RESIDUAL_H
#define NP_REGRESSION_ALLLARGE_RESIDUAL_H

/* Accurate uncertainty moments for the incumbent all-large polynomial map.
 * These bounded expansions retain the supplied double basis/inverse exactly;
 * they do not refactor, replace or regularize the point-estimation solve.
 * The ordinary representation reuses the paired-contrast expansion arithmetic.
 * Only an explicit exponent-representability failure admits the cold tagged
 * representation. Capacity/nonfinite failures remain failures in both modes.
 */
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "regression_contrast.h"

enum { NP_ALLLARGE_OK=0, NP_ALLLARGE_INVALID=1, NP_ALLLARGE_EXPONENT=2 };

static inline int np_alllarge_error(int a,int b)
{
  return a==NP_ALLLARGE_INVALID||b==NP_ALLLARGE_INVALID ? NP_ALLLARGE_INVALID :
    (a||b ? NP_ALLLARGE_EXPONENT : NP_ALLLARGE_OK);
}

/* Same grow-expansion operation as np_contrast_exact_add, with a distinct
 * exponent failure so capacity exhaustion cannot silently select a retry. */
static inline void np_alllarge_narrow_add_double(NPContrastExactSum *a,double value)
{
  if(a->invalid)return;
  if(!isfinite(value)){a->invalid=NP_ALLLARGE_INVALID;return;}
  double out[NP_CONTRAST_EXACT_CAPACITY];int count=0;
  for(int i=0;i<a->size;++i){
    const NPContrastNumber part=np_contrast_two_sum(value,a->component[i]);
    if(!np_contrast_finite(part)){a->invalid=NP_ALLLARGE_EXPONENT;return;}
    if(part.lo!=0.0){
      if(count==NP_CONTRAST_EXACT_CAPACITY){a->invalid=NP_ALLLARGE_INVALID;return;}
      out[count++]=part.lo;
    }
    value=part.hi;
  }
  if(value!=0.0){
    if(count==NP_CONTRAST_EXACT_CAPACITY){a->invalid=NP_ALLLARGE_INVALID;return;}
    out[count++]=value;
  }
  memcpy(a->component,out,(size_t)count*sizeof(*out));a->size=count;
}

static inline void np_alllarge_narrow_compress(NPContrastExactSum *a)
{
  if(a->invalid||a->size<2)return;
  NPContrastExactSum out={0};
  for(int i=a->size-1;i>=0;--i)np_alllarge_narrow_add_double(&out,a->component[i]);
  *a=out;
}

static inline void np_alllarge_narrow_product_add(NPContrastExactSum *a,double x,double y)
{
  if(a->invalid)return;
  if(!isfinite(x)||!isfinite(y)){a->invalid=NP_ALLLARGE_INVALID;return;}
  if(x==0.0||y==0.0)return;
  NPContrastNumber product;
  if(fabs(x)>=0x1p-450&&fabs(x)<=0x1p450&&
     fabs(y)>=0x1p-450&&fabs(y)<=0x1p450){
    product=np_contrast_two_product(x,y);
  }else{
    int xe,ye;
    const double xm=frexp(x,&xe),ym=frexp(y,&ye);
    const NPContrastNumber m=np_contrast_two_product(xm,ym);
    product.hi=ldexp(m.hi,xe+ye);product.lo=ldexp(m.lo,xe+ye);
    if(!isfinite(product.hi)||!isfinite(product.lo)||
       ldexp(product.hi,-xe-ye)!=m.hi||ldexp(product.lo,-xe-ye)!=m.lo){
      a->invalid=NP_ALLLARGE_EXPONENT;return;
    }
  }
  np_alllarge_narrow_add_double(a,product.lo);
  np_alllarge_narrow_add_double(a,product.hi);
}

static inline NPContrastExactSum np_alllarge_narrow_scale(const NPContrastExactSum *a,double x)
{
  NPContrastExactSum out={0};out.invalid=a->invalid;
  for(int i=0;i<a->size&&!out.invalid;++i)
    np_alllarge_narrow_product_add(&out,a->component[i],x);
  return out;
}

static inline void np_alllarge_narrow_add(NPContrastExactSum *a,const NPContrastExactSum *b,double sign)
{
  a->invalid=np_alllarge_error(a->invalid,b->invalid);
  for(int i=0;i<b->size&&!a->invalid;++i){
    if(sign==1.0||sign==-1.0)
      np_alllarge_narrow_add_double(a,sign*b->component[i]);
    else
      np_alllarge_narrow_product_add(a,b->component[i],sign);
  }
}

static inline NPContrastExactSum np_alllarge_narrow_multiply(const NPContrastExactSum *a,const NPContrastExactSum *b)
{
  NPContrastExactSum out={0};out.invalid=np_alllarge_error(a->invalid,b->invalid);
  for(int i=0;i<a->size&&!out.invalid;++i)for(int j=0;j<b->size;++j)
    np_alllarge_narrow_product_add(&out,a->component[i],b->component[j]);
  return out;
}

static inline int np_alllarge_narrow_direction_sign(const NPContrastExactSum *a)
{
  return a->size==0?0:(a->component[a->size-1]>0.0?1:-1);
}

static inline double np_alllarge_narrow_mantissa(const NPContrastExactSum *a,int *exponent)
{
  if(a->size==0){*exponent=0;return 0.0;}
  (void)frexp(a->component[a->size-1],exponent);
  double value=0.0;
  for(int i=0;i<a->size;++i)value+=ldexp(a->component[i],-*exponent);
  return value;
}

/* Cold components carry exponents separately. Multiplication always acts on
 * finite normalized mantissas; finite input-double expression exponents fit
 * in int. A component farther than 54 bits cannot overlap the other operand,
 * so retaining it separately is exact, not an underflow cutoff. */
typedef struct {double mantissa;int exponent;} NPAllLargeWideNumber;
typedef struct {
  NPAllLargeWideNumber component[NP_CONTRAST_EXACT_CAPACITY];
  int size,invalid;
} NPAllLargeWideExpansion;

static inline NPAllLargeWideNumber np_alllarge_wide_number(double value,int exponent)
{
  NPAllLargeWideNumber out={0};
  if(value!=0.0){int shift;out.mantissa=frexp(value,&shift);out.exponent=exponent+shift;}
  return out;
}

static inline void np_alllarge_wide_two_sum(NPAllLargeWideNumber a,NPAllLargeWideNumber b,
  NPAllLargeWideNumber *sum,NPAllLargeWideNumber *error)
{
  *error=(NPAllLargeWideNumber){0};
  if(a.mantissa==0.0){*sum=b;return;}
  if(b.mantissa==0.0){*sum=a;return;}
  if(a.exponent<b.exponent){const NPAllLargeWideNumber swap=a;a=b;b=swap;}
  const int gap=a.exponent-b.exponent;
  if(gap>54){*sum=a;*error=b;return;}
  const NPContrastNumber joined=np_contrast_two_sum(a.mantissa,ldexp(b.mantissa,-gap));
  *sum=np_alllarge_wide_number(joined.hi,a.exponent);
  *error=np_alllarge_wide_number(joined.lo,a.exponent);
}

static inline void np_alllarge_wide_add_number(NPAllLargeWideExpansion *a,NPAllLargeWideNumber value)
{
  if(a->invalid||value.mantissa==0.0)return;
  if(!isfinite(value.mantissa)){a->invalid=NP_ALLLARGE_INVALID;return;}
  NPAllLargeWideNumber out[NP_CONTRAST_EXACT_CAPACITY];int count=0;
  for(int i=0;i<a->size;++i){
    NPAllLargeWideNumber sum,error;
    np_alllarge_wide_two_sum(value,a->component[i],&sum,&error);
    if(error.mantissa!=0.0){
      if(count==NP_CONTRAST_EXACT_CAPACITY){a->invalid=NP_ALLLARGE_INVALID;return;}
      out[count++]=error;
    }
    value=sum;
  }
  if(value.mantissa!=0.0){
    if(count==NP_CONTRAST_EXACT_CAPACITY){a->invalid=NP_ALLLARGE_INVALID;return;}
    out[count++]=value;
  }
  memcpy(a->component,out,(size_t)count*sizeof(*out));a->size=count;
}

static inline void np_alllarge_wide_add_double(NPAllLargeWideExpansion *a,double value)
{
  if(!isfinite(value)){a->invalid=NP_ALLLARGE_INVALID;return;}
  np_alllarge_wide_add_number(a,np_alllarge_wide_number(value,0));
}

static inline void np_alllarge_wide_compress(NPAllLargeWideExpansion *a)
{
  if(a->invalid||a->size<2)return;
  NPAllLargeWideExpansion out={0};
  for(int i=a->size-1;i>=0;--i)np_alllarge_wide_add_number(&out,a->component[i]);
  *a=out;
}

static inline void np_alllarge_wide_product_numbers(NPAllLargeWideExpansion *a,
  NPAllLargeWideNumber x,NPAllLargeWideNumber y)
{
  if(x.mantissa==0.0||y.mantissa==0.0)return;
  const NPContrastNumber product=np_contrast_two_product(x.mantissa,y.mantissa);
  np_alllarge_wide_add_number(a,np_alllarge_wide_number(product.lo,x.exponent+y.exponent));
  np_alllarge_wide_add_number(a,np_alllarge_wide_number(product.hi,x.exponent+y.exponent));
}

static inline void np_alllarge_wide_product_add(NPAllLargeWideExpansion *a,double x,double y)
{
  if(!isfinite(x)||!isfinite(y)){a->invalid=NP_ALLLARGE_INVALID;return;}
  np_alllarge_wide_product_numbers(a,np_alllarge_wide_number(x,0),np_alllarge_wide_number(y,0));
}

static inline NPAllLargeWideExpansion np_alllarge_wide_scale(const NPAllLargeWideExpansion *a,double x)
{
  NPAllLargeWideExpansion out={0};out.invalid=a->invalid||!isfinite(x);
  const NPAllLargeWideNumber factor=np_alllarge_wide_number(x,0);
  for(int i=0;i<a->size&&!out.invalid;++i)
    np_alllarge_wide_product_numbers(&out,a->component[i],factor);
  return out;
}

static inline void np_alllarge_wide_add(NPAllLargeWideExpansion *a,const NPAllLargeWideExpansion *b,double sign)
{
  a->invalid=np_alllarge_error(a->invalid,b->invalid);
  for(int i=0;i<b->size&&!a->invalid;++i){
    if(sign==1.0||sign==-1.0){
      NPAllLargeWideNumber value=b->component[i];value.mantissa*=sign;
      np_alllarge_wide_add_number(a,value);
    }else{
      np_alllarge_wide_product_numbers(a,b->component[i],np_alllarge_wide_number(sign,0));
    }
  }
}

static inline NPAllLargeWideExpansion np_alllarge_wide_multiply(const NPAllLargeWideExpansion *a,const NPAllLargeWideExpansion *b)
{
  NPAllLargeWideExpansion out={0};out.invalid=np_alllarge_error(a->invalid,b->invalid);
  for(int i=0;i<a->size&&!out.invalid;++i)for(int j=0;j<b->size;++j)
    np_alllarge_wide_product_numbers(&out,a->component[i],b->component[j]);
  return out;
}

static inline int np_alllarge_wide_direction_sign(const NPAllLargeWideExpansion *a)
{
  return a->size==0?0:(a->component[a->size-1].mantissa>0.0?1:-1);
}

static inline double np_alllarge_wide_mantissa(const NPAllLargeWideExpansion *a,int *exponent)
{
  if(a->size==0){*exponent=0;return 0.0;}
  *exponent=a->component[a->size-1].exponent;
  double value=0.0;
  for(int i=0;i<a->size;++i)
    value+=ldexp(a->component[i].mantissa,a->component[i].exponent-*exponent);
  return value;
}

#define NP_AL_JOIN_IMPL(a,b) a##b
#define NP_AL_JOIN(a,b) NP_AL_JOIN_IMPL(a,b)
#define NP_AL_NAME(name) NP_AL_JOIN(np_alllarge_narrow_,name)
#define NP_AL_EXPANSION NPContrastExactSum
#include "regression_alllarge_moments.h"
#undef NP_AL_NAME
#undef NP_AL_EXPANSION
#define NP_AL_NAME(name) NP_AL_JOIN(np_alllarge_wide_,name)
#define NP_AL_EXPANSION NPAllLargeWideExpansion
#include "regression_alllarge_moments.h"
#undef NP_AL_NAME
#undef NP_AL_EXPANSION
#undef NP_AL_JOIN
#undef NP_AL_JOIN_IMPL

typedef struct {
  np_alllarge_narrow_moments narrow;
  np_alllarge_wide_moments wide;
} NPAllLargeResidualWorkspace;

static inline void np_alllarge_residual_workspace_clear(NPAllLargeResidualWorkspace *workspace)
{
  if(workspace!=NULL){
    np_alllarge_narrow_clear(&workspace->narrow);
    np_alllarge_wide_clear(&workspace->wide);
  }
}
#endif
