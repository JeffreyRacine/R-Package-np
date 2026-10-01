#ifndef NP_CONDITIONAL_DELETED_DOT_H
#define NP_CONDITIONAL_DELETED_DOT_H
/* Same twofold sum/np_cqr_product, four independent dependency chains. No fast-math. */
static double np_cqr_dot4(const double *x,int sx,const double *y,int sy,int n){
  NPConditionalQRAcc a={0,0},b={0,0},c={0,0},d={0,0};int i=0;
  for(;i+3<n;i+=4){
    np_cqr_product(&a,x[(size_t)i*sx],y[(size_t)i*sy]);
    np_cqr_product(&b,x[(size_t)(i+1)*sx],y[(size_t)(i+1)*sy]);
    np_cqr_product(&c,x[(size_t)(i+2)*sx],y[(size_t)(i+2)*sy]);
    np_cqr_product(&d,x[(size_t)(i+3)*sx],y[(size_t)(i+3)*sy]);
  }
  for(;i<n;i++)np_cqr_product(&a,x[(size_t)i*sx],y[(size_t)i*sy]);
  np_cqr_add(&a,b.hi);np_cqr_add(&a,b.lo);np_cqr_add(&a,c.hi);np_cqr_add(&a,c.lo);np_cqr_add(&a,d.hi);np_cqr_add(&a,d.lo);
  return a.hi+a.lo;
}

#endif
