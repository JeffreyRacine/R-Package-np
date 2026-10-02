/* Conditional GNN coefficients: twofold triangular arithmetic and local polynomial
 * evaluation. Factorization and rank decisions are unchanged. */
#include "conditional_deleted_qr.h"
static NPConditionalQRAcc cqr_dd_mul(NPConditionalQRAcc a,NPConditionalQRAcc b){
  NPConditionalQRAcc out={0,0};
  np_cqr_product(&out,a.hi,b.hi);np_cqr_product(&out,a.hi,b.lo);
  np_cqr_product(&out,a.lo,b.hi);np_cqr_product(&out,a.lo,b.lo);return out;
}
static NPConditionalQRAcc cqr_dd_div(NPConditionalQRAcc a,double b){
  const double q=a.hi/b;
  np_cqr_product(&a,-q,b);
  NPConditionalQRAcc out={q,0};np_cqr_add(&out,(a.hi+a.lo)/b);return out;
}
static int cqr_compensated_row(const NPConditionalQRFast *f,int n,int k,
    const int *terms,const double *x,const double *w,int pos,
    NPConditionalQRAcc *v,NPConditionalQRAcc *coef,double *row){
  double maxw=0,radius=0;
  for(int i=0;i<n;++i)if(i!=pos && w[i]>0){maxw=fmax(maxw,w[i]);radius=fmax(radius,fabs(x[i]-x[pos]));}
  if(!(maxw>0 && radius>0))return 1;
  for(int i=0;i<k;++i){
    NPConditionalQRAcc value={f->pivot[i]==1?1.:0.,0};
    value=cqr_dd_div(value,f->scale[0]);
    for(int j=0;j<i;++j){
      np_cqr_product(&value,-f->a[j+(size_t)f->capacity*i],v[j].hi);
      np_cqr_product(&value,-f->a[j+(size_t)f->capacity*i],v[j].lo);
    }
    v[i]=cqr_dd_div(value,f->a[i+(size_t)f->capacity*i]);
  }
  for(int i=k-1;i>=0;--i){
    NPConditionalQRAcc value=v[i];
    for(int j=i+1;j<k;++j){
      np_cqr_product(&value,-f->a[i+(size_t)f->capacity*j],v[j].hi);
      np_cqr_product(&value,-f->a[i+(size_t)f->capacity*j],v[j].lo);
    }
    v[i]=cqr_dd_div(value,f->a[i+(size_t)f->capacity*i]);
  }
  for(int i=0;i<k;++i){int p=f->pivot[i]-1;coef[p]=cqr_dd_div(cqr_dd_div(v[i],f->scale[p]),maxw);}
  for(int i=0;i<n;++i){
    if(i==pos || w[i]==0){row[i]=0;continue;}
    NPConditionalQRAcc u={x[i],0},sum={0,0};np_cqr_add(&u,-x[pos]);u=cqr_dd_div(u,radius);
    for(int t=0;t<k;++t){
      NPConditionalQRAcc power={1,0};
      for(int p=0;p<terms[t];++p)power=cqr_dd_mul(power,u);
      NPConditionalQRAcc term=cqr_dd_mul(coef[t],power);
      np_cqr_add(&sum,term.hi);np_cqr_add(&sum,term.lo);
    }
    NPConditionalQRAcc weight={w[i],0};sum=cqr_dd_mul(sum,weight);row[i]=sum.hi+sum.lo;
    if(!R_FINITE(row[i]))return 1;
  }
  return 0;
}
