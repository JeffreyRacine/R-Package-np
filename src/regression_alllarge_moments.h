/* Intentional arithmetic template: included once for ordinary bounded exact
 * expansions and once for the cold exponent-tagged representation. Only the
 * arithmetic representation changes; the retained-map identity is identical.
 * G=B'B, beta=A B'Y, g=A B'1 + D00 A[,0], M=A G A'. For row b, t=b'g:
 * e=tY-b'beta; q=t^2-2t(b'Ab)+b'Mb. This equals the admitted off-diagonal
 * completion with residual diagonal D00*(b'A)[0]+sum_off(S), including ridge.
 */
typedef struct {
  int p;
  NP_AL_EXPANSION *allocation,*gram,*meat,*coefficient,*total,*cross,*sums,*temporary;
  size_t count;
} NP_AL_NAME(moments);

static inline void NP_AL_NAME(clear)(NP_AL_NAME(moments) *context)
{
  free(context->allocation);memset(context,0,sizeof(*context));
}

static inline int NP_AL_NAME(prepare)(NP_AL_NAME(moments) *context,int n,int p,
  double **basis,const double *response,const double *inverse,
  NPContrastNumber intercept_ridge,void (*activity)(void))
{
  if(p<=0||n<=0||(size_t)p>SIZE_MAX/(size_t)p)return NP_ALLLARGE_INVALID;
  const size_t square=(size_t)p*p;
  if(square>SIZE_MAX/3||(size_t)p>(SIZE_MAX-3*square)/4)return NP_ALLLARGE_INVALID;
  const size_t count=3*square+4*(size_t)p;
  if(count>SIZE_MAX/sizeof(NP_AL_EXPANSION))return NP_ALLLARGE_INVALID;
  context->allocation=(NP_AL_EXPANSION *)calloc(count,sizeof(NP_AL_EXPANSION));
  if(context->allocation==NULL)return NP_ALLLARGE_INVALID;
  context->count=count;context->p=p;
  context->gram=context->allocation;context->meat=context->gram+square;
  context->temporary=context->meat+square;
  context->coefficient=context->temporary+square;
  context->total=context->coefficient+p;context->cross=context->total+p;
  context->sums=context->cross+p;
  size_t work=0;
  for(int i=0;i<n;++i)for(int a=0;a<p;++a){
    NP_AL_NAME(add_double)(&context->sums[a],basis[a][i]);
    NP_AL_NAME(product_add)(&context->cross[a],basis[a][i],response[i]);
    for(int b=0;b<=a;++b){
      NP_AL_NAME(product_add)(&context->gram[(size_t)a*p+b],basis[a][i],basis[b][i]);
      if(activity!=NULL&&((++work)&4095)==0)activity();
    }
  }
  for(int a=0;a<p;++a)for(int b=0;b<a;++b)
    context->gram[(size_t)b*p+a]=context->gram[(size_t)a*p+b];
  for(int a=0;a<p;++a){
    for(int k=0;k<p;++k){
      NP_AL_EXPANSION value=NP_AL_NAME(scale)(&context->cross[k],inverse[a+(size_t)k*p]);
      NP_AL_NAME(add)(&context->coefficient[a],&value,1.0);
      value=NP_AL_NAME(scale)(&context->sums[k],inverse[a+(size_t)k*p]);
      NP_AL_NAME(add)(&context->total[a],&value,1.0);
      for(int b=0;b<p;++b){
        value=NP_AL_NAME(scale)(&context->gram[(size_t)k*p+b],inverse[a+(size_t)k*p]);
        NP_AL_NAME(add)(&context->temporary[(size_t)a*p+b],&value,1.0);
        if(activity!=NULL&&((++work)&4095)==0)activity();
      }
    }
    NP_AL_NAME(product_add)(&context->total[a],intercept_ridge.hi,inverse[a]);
    NP_AL_NAME(product_add)(&context->total[a],intercept_ridge.lo,inverse[a]);
  }
  for(int a=0;a<p;++a)for(int b=0;b<=a;++b){
    for(int k=0;k<p;++k){
      NP_AL_EXPANSION value=NP_AL_NAME(scale)(&context->temporary[(size_t)a*p+k],inverse[b+(size_t)k*p]);
      NP_AL_NAME(add)(&context->meat[(size_t)a*p+b],&value,1.0);
      if(activity!=NULL&&((++work)&4095)==0)activity();
    }
    context->meat[(size_t)b*p+a]=context->meat[(size_t)a*p+b];
  }
  int status=NP_ALLLARGE_OK;
  for(size_t i=0;i<count;++i){
    status=np_alllarge_error(status,context->allocation[i].invalid);
    if(activity!=NULL&&((++work)&4095)==0)activity();
  }
  if(status!=NP_ALLLARGE_OK)return status;
  for(size_t i=0;i<square;++i){
    NP_AL_NAME(compress)(&context->gram[i]);NP_AL_NAME(compress)(&context->meat[i]);
  }
  for(int a=0;a<p;++a){
    NP_AL_NAME(compress)(&context->coefficient[a]);NP_AL_NAME(compress)(&context->total[a]);
  }
  for(size_t i=0;i<count;++i)
    status=np_alllarge_error(status,context->allocation[i].invalid);
  if(activity!=NULL)activity();
  return status;
}

static inline int NP_AL_NAME(row)(const NP_AL_NAME(moments) *context,double **basis,
  int self,double response,const double *inverse,double *normalized,
  int *information,void (*activity)(void))
{
  const int p=context->p;
  NP_AL_EXPANSION diagonal={0},quadratic={0},total={0},fitted={0};
  size_t work=0;
  for(int a=0;a<p;++a){
    NP_AL_EXPANSION value=NP_AL_NAME(scale)(&context->coefficient[a],basis[a][self]);
    NP_AL_NAME(add)(&fitted,&value,1.0);
    value=NP_AL_NAME(scale)(&context->total[a],basis[a][self]);
    NP_AL_NAME(add)(&total,&value,1.0);
    for(int k=0;k<p;++k){
      NP_AL_EXPANSION part={0};
      NP_AL_NAME(product_add)(&part,basis[a][self],inverse[a+(size_t)k*p]);
      value=NP_AL_NAME(scale)(&part,basis[k][self]);NP_AL_NAME(add)(&diagonal,&value,1.0);
      value=NP_AL_NAME(scale)(&context->meat[(size_t)a*p+k],basis[a][self]);
      part=NP_AL_NAME(scale)(&value,basis[k][self]);NP_AL_NAME(add)(&quadratic,&part,1.0);
      if(activity!=NULL&&((++work)&4095)==0)activity();
    }
  }
  NP_AL_EXPANSION variance=quadratic;
  NP_AL_EXPANSION value=NP_AL_NAME(multiply)(&total,&total);
  NP_AL_NAME(add)(&variance,&value,1.0);
  value=NP_AL_NAME(multiply)(&total,&diagonal);NP_AL_NAME(add)(&variance,&value,-2.0);
  NP_AL_EXPANSION numerator=NP_AL_NAME(scale)(&total,response);
  NP_AL_NAME(add)(&numerator,&fitted,-1.0);
  const int status=np_alllarge_error(variance.invalid,numerator.invalid);
  if(status!=NP_ALLLARGE_OK)return status;
  const int sign=NP_AL_NAME(direction_sign)(&variance);
  if(sign<0)return NP_ALLLARGE_INVALID;
  *normalized=0.0;
  if(sign==0){
    if(NP_AL_NAME(direction_sign)(&numerator)!=0)return NP_ALLLARGE_INVALID;
    *information=NP_RESIDUAL_UNIDENTIFIED;
    return NP_ALLLARGE_OK;
  }
  /* Never convert q itself to double: a positive norm may be below the
   * subnormal range while its normalized residual is ordinary finite. */
  int qe,ee;
  double qm=NP_AL_NAME(mantissa)(&variance,&qe);
  const double em=NP_AL_NAME(mantissa)(&numerator,&ee);
  if(!(qm>0.0)||!isfinite(qm)||!isfinite(em))return NP_ALLLARGE_INVALID;
  if(qe%2!=0){qm*=2.0;--qe;}
  *normalized=ldexp(em/sqrt(qm),ee-qe/2);
  if(!isfinite(*normalized))return NP_ALLLARGE_INVALID;
  *information=NP_RESIDUAL_IDENTIFIED;
  return NP_ALLLARGE_OK;
}
