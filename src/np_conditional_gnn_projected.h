/* Bounded conditional scalar-response GNN contraction, extracted from the
 * preserved qualified C167 computation. Included after canonical row providers.
 * The unconditional owner and compact-X prefix implementation are untouched. */
typedef int (*NPGNNIntegralFilledRow)(void *,int,double *,double *);
typedef struct {
  int n,folded,has_categories,nprofiles,first_fold,end_fold,pair_rank,pair_ranks,kernel;
  double logG[17];
  double **data,*values,*bounds,*weight_row,*left,*right,*overlap,*products;
  double *profile_moments,*profile_left,*profile_right,*pair_rows;
  int *profile_id,*profile_rep;
  size_t budget;
  NPGNNIntegralGeometry *geometry;
  NPGNNIntegralWeightRow row;
  NPGNNIntegralFilledRow filled_row;
  NPGNNIntegralCategoryOverlap category_overlap;
  void *context;
  double capped;
  int *projected_order,*projected_start;
  double *projected_factor,*projected_diagonal,*projected_residual;
  double *projected_recursive,*projected_estimates;
  double *compact_cuts;
  int compact_pair_geometry;
  int *projected_selected;
  int projected_rows,projected_intervals,projected_shared;
  size_t projected_bytes;
} NPGNNConditionalProjectionOwner;

static int np_gnn_integral_reciprocal_interval(double qlo, double qhi, double anchor,
                                               double *lo, double *hi, double *width)
{
  const double u0 = qlo == anchor ? R_PosInf : 1.0/(qlo-anchor);
  const double u1 = qhi == anchor ? R_NegInf : 1.0/(qhi-anchor);
  *lo = fmin(u0,u1); *hi = fmax(u0,u1); *width = *hi-*lo;
  if(*lo < *hi) return 0;
  if(*lo != *hi || !R_FINITE(*lo) || !(qlo < qhi)) return 1;
  const long double a = (long double)qlo-anchor, b = (long double)qhi-anchor;
  if(a == 0 || b == 0) return 1;
  *width = (double)fabsl(((long double)qhi-qlo)/a/b);
  return !R_FINITE(*width) || *width < 0;
}

/* contract bounded influence rows before squaring. */
static int np_gnn_projected_storage(NPGNNConditionalProjectionOwner *c)
{
  const size_t n=(size_t)c->n;
  const size_t compact_bytes=c->kernel>=4 ? (2*n+2)*sizeof(double) : 0;
  /* Shared-basis storage is linear in n with at most 512 seed columns. */
  if(c->folded && !c->has_categories && c->filled_row!=NULL &&
     c->first_fold==0 && c->end_fold==c->n) {
    size_t bytes=sizeof(*c)+sizeof(NPGNNIntegralGeometry)+compact_bytes,cells;
    int invalid=n>(size_t)INT_MAX/192 ||
       np_gnn_integral_size_add(&bytes,n,4*sizeof(NPGNNIntegralInterval)+
         4*sizeof(NPGNNIntegralPiece)+12*sizeof(double)+4*sizeof(int)) ||
       !np_size_mul_checked(n,4*512+3*192,&cells) ||
       np_gnn_integral_size_add(&bytes,cells,sizeof(double)) ||
       np_gnn_integral_size_add(&bytes,192*512,sizeof(double)) ||
       !np_size_mul_checked(n,3*21+2*16,&cells) ||
       np_gnn_integral_size_add(&bytes,cells,sizeof(double)) ||
       np_gnn_integral_size_add(&bytes,n+512,sizeof(int)) ||
       np_gnn_integral_size_add(&bytes,512+1024,sizeof(double));
    if(!invalid && bytes<=c->budget) {
      c->projected_rows=(int)n;c->projected_intervals=8;c->projected_bytes=bytes;
      c->projected_shared=1;

      return 1;
    }

  }
  const size_t groups=c->has_categories?(size_t)c->nprofiles:1;
  int rows=MIN(512,c->folded?c->n:1),intervals=8;
  size_t fixed=compact_bytes,cells;
  if(groups<1 ||
     np_gnn_integral_size_add(&fixed,1,sizeof(*c)+sizeof(NPGNNIntegralGeometry)) ||
     np_gnn_integral_size_add(&fixed,n,4*sizeof(NPGNNIntegralInterval)+
       4*sizeof(NPGNNIntegralPiece)+sizeof(double)+sizeof(int)) ||
     np_gnn_integral_size_add(&fixed,n,3*sizeof(int)+4*sizeof(double)) ||
     np_gnn_integral_size_add(&fixed,groups+1,sizeof(int)) ||
     !np_size_mul_checked(groups,groups,&cells) ||
     np_gnn_integral_size_add(&fixed,cells,sizeof(double)) ||
     /* Account for maximum-depth recursion and initial record scratch. */
     np_gnn_integral_size_add(&fixed,(size_t)512*(3*21+2*16),sizeof(double)) ||
     np_gnn_integral_size_add(&fixed,512,sizeof(int)) ||
     np_gnn_integral_size_add(&fixed,1024,sizeof(double)))return 0;
  for(;;) {
    const size_t nodes=(size_t)intervals*24;
    size_t bytes=fixed;
    int invalid=n>(size_t)INT_MAX/nodes || groups>(size_t)INT_MAX/nodes ||
      !np_size_mul_checked(n,nodes,&cells) ||
      np_gnn_integral_size_add(&bytes,cells,2*sizeof(double)) ||
      !np_size_mul_checked(n,(size_t)rows,&cells) ||
      np_gnn_integral_size_add(&bytes,cells,sizeof(double)) ||
      !np_size_mul_checked(nodes,(size_t)rows,&cells) ||
      !np_size_mul_checked(cells,groups,&cells) ||
      np_gnn_integral_size_add(&bytes,cells,sizeof(double));
    if(!invalid && bytes<=c->budget) {
      c->projected_rows=rows;c->projected_intervals=intervals;
      c->projected_bytes=bytes;
      return 1;
    }
    if(rows>1)rows=(rows+1)/2;
    else if(intervals>1)intervals=(intervals+1)/2;
    else return 0;
  }
}

typedef struct {
  NPGNNConditionalProjectionOwner *owner;
  int rows,groups,rank,first,*selected;
  double anchor,norm_max,category_max,absolute_scale,compression_bound;
} NPGNNProjectedIntegral;

/* Convex-hull bounds of the canonical polynomials in Bernstein form on
 * z^2 in [0,5]. Maximal absolute coefficients before the kernel multiplier
 * are 1, 75, 175/48 and 3675/512, respectively. The slack covers the small
 * fixed polynomial evaluation, not quadrature error or a changed tolerance. */
static double np_gnn_projected_compact_height(int kernel)
{
  const double base=.33541019662496845446;
  const double height=kernel==8 ? .5 : kernel==4 ? base :
    kernel==5 ? .008385254916*75 : kernel==6 ? base*(175.0/48) :
    base*(3675.0/512);
  return nextafter(height*(1+256*DBL_EPSILON),R_PosInf);
}

/* Fill exactly one bounded slab in the same Y-sorted fold order. The
 * ordinary deleted row is never overwritten outside this private slab. */
static int np_gnn_projected_fill_slab(NPGNNProjectedIntegral *q,
  int first,int width)
{
  NPGNNConditionalProjectionOwner *c=q->owner;
  const int n=c->n;
  for(int f=0;f<width;++f) {
    if((f & 15)==0)np_progress_bandwidth_loop_step();
    const int fold=c->geometry->order[first+f];
    double diagonal=R_NaN;
    if(c->filled_row(c->context,fold,c->weight_row,&diagonal) ||
       !R_FINITE(diagonal) || c->weight_row[fold]!=0.0)return 1;
    c->projected_diagonal[first+f]=diagonal;
    c->weight_row[fold]=diagonal;
    double norm=0;
    for(int i=0;i<n;++i) {
      const double value=c->weight_row[c->projected_order[i]];
      if(!R_FINITE(value))return 1;
      c->right[(size_t)f*n+i]=value;norm+=fabs(value);
    }
    q->norm_max=fmax(q->norm_max,norm);
  }
  return 0;
}

/* Compensated residual of stored Q/R. FMA retains the
 * multiplication residual; TwoSum retains the addition residual. The bound
 * includes accumulation of the small corrections and underflow slack. */
static double np_gnn_projected_residual(const double *q,int stride,
  const double *r,int rank,double original)
{
  double sum=0,correction=0,correction_mass=0;
  for(int j=0;j<rank;++j) {
    const double product=q[(size_t)j*stride]*r[j];
    const double product_error=fma(q[(size_t)j*stride],r[j],-product);
    const double next=sum+product,back=next-sum;
    const double addition_error=(sum-(next-back))+(product-back);
    correction+=product_error;
    correction+=addition_error;
    correction_mass+=fabs(product_error)+fabs(addition_error);
    sum=next;
  }
  const double difference=original-sum,residual=difference-correction;
  const double gamma=(4*rank+8)*DBL_EPSILON;
  return nextafter(fabs(residual)+gamma*correction_mass+
    4*DBL_EPSILON*(fabs(difference)+fabs(residual))+(rank+1)*DBL_MIN,R_PosInf);
}

/* Certified QR compression of the filled influence plane. The
 * removed diagonal is exact and separate; no deletion is approximated. */
static int np_gnn_projected_factor(NPGNNProjectedIntegral *q)
{
  NPGNNConditionalProjectionOwner *c=q->owner;
  const int n=c->n,b=MIN(512,q->rows),p=MIN(n,b),folds=q->rows;
  if(q->groups!=1 || c->has_categories || !c->folded)return 1;
  if(np_gnn_projected_fill_slab(q,0,b))return 1;
  memcpy(c->projected_factor,c->right,(size_t)n*b*sizeof(double));
  int pivots[512]={0},info=0,lwork=-1;
  double tau[512],query;
  F77_CALL(dgeqp3)(&n,&b,c->projected_factor,&n,pivots,tau,&query,&lwork,&info);
  if(info || !R_FINITE(query) || query<1 || query>INT_MAX)return 1;
  lwork=(int)ceil(query);

  if((size_t)lwork>(c->budget-c->projected_bytes)/sizeof(double)) {
    return 2;
  }
  double *work=(double *)malloc((size_t)lwork*sizeof(double));
  if(!work)return 1;
  np_progress_bandwidth_loop_step();
  F77_CALL(dgeqp3)(&n,&b,c->projected_factor,&n,pivots,tau,work,&lwork,&info);
  np_progress_bandwidth_loop_step();
  free(work);if(info)return 1;
  int rank=p;
  long double tail=0;
  for(int row=p-1;row>0;--row) {
    for(int col=row;col<b;++col) {
      const double value=c->projected_factor[(size_t)col*n+row];
      tail+=(long double)value*value;
    }
    if(tail<1e-30L)rank=row;else break;
  }
  rank=MIN(p,rank+4);
  c->profile_left=(double *)calloc((size_t)rank*folds,sizeof(double));
  if(!c->profile_left)return 1;
  lwork=-1;
  F77_CALL(dorgqr)(&n,&rank,&rank,c->projected_factor,&n,tau,&query,&lwork,&info);
  if(info || !R_FINITE(query) || query<1 || query>INT_MAX)return 1;
  lwork=(int)ceil(query);
  if((size_t)lwork>(c->budget-c->projected_bytes)/sizeof(double)) {
    return 2;
  }
  work=(double *)malloc((size_t)lwork*sizeof(double));
  if(!work)return 1;
  np_progress_bandwidth_loop_step();
  F77_CALL(dorgqr)(&n,&rank,&rank,c->projected_factor,&n,tau,work,&lwork,&info);
  np_progress_bandwidth_loop_step();
  free(work);if(info)return 1;
  /* The representation certificate uses a conservative
   * whole-domain envelope, independent of quadrature samples. */
  long double reciprocal_width=0,max_bound=0;
  double minimum_bandwidth=R_PosInf;
  for(int z=0;z<c->geometry->count;++z) {
    if((z & 4095)==0)np_progress_bandwidth_loop_step();
    const NPGNNIntegralInterval t=c->geometry->intervals[z];
    double lo,hi,width0,width1;
    if(np_gnn_integral_reciprocal_interval(t.lo,t.hi,t.primary_anchor,&lo,&hi,&width0) ||
       np_gnn_integral_reciprocal_interval(t.lo,t.hi,t.successor_anchor,&lo,&hi,&width1))return 1;
    reciprocal_width+=fmax(width0,width1);
    minimum_bandwidth=fmin(minimum_bandwidth,c->geometry->scale*
      fmin(fmin(fabs(t.lo-t.primary_anchor),fabs(t.hi-t.primary_anchor)),
           fmin(fabs(t.lo-t.successor_anchor),fabs(t.hi-t.successor_anchor))));
  }
  const long double height=(c->kernel==0 ? allck[0](0) : c->kernel<4 ?
    exp(c->logG[0]) : np_gnn_projected_compact_height(c->kernel))/c->geometry->scale;
  /* For any h(q)>=hmin, K_h(q-y)^2 is bounded by its pointwise
   * supremum over h>=hmin. Integrating that envelope over the whole line
   * gives [integral(-1,1) phi(t)^2 dt + 2 phi(1)^2]/hmin. */
  const double envelope_constant=allck[0](0)/sqrt(2.0)*
    (1-2*pnorm5(-sqrt(2.0),0,1,1,0))+2*allck[0](1)*allck[0](1);
  const long double integral_bound=c->kernel==0 ?
    fminl(reciprocal_width*height*height,
      nextafter(envelope_constant/minimum_bandwidth,R_PosInf)) :
    reciprocal_width*height*height;
  const double one=1,zero=0,minus=-1;
  for(int first=0;first<folds;first+=512) {
    np_progress_bandwidth_loop_step();
    const int width=MIN(512,folds-first);
    if(first && np_gnn_projected_fill_slab(q,first,width))return 1;
    double *coeff=c->profile_left+(size_t)first*rank;
    F77_CALL(dgemm)("T","N",&rank,&width,&n,&one,c->projected_factor,&n,
      c->right,&n,&zero,coeff,&rank FCONE FCONE);
    memcpy(c->projected_residual,c->right,(size_t)n*width*sizeof(double));
    F77_CALL(dgemm)("N","N",&n,&width,&rank,&minus,c->projected_factor,&n,
      coeff,&rank,&one,c->projected_residual,&n FCONE FCONE);
    F77_CALL(dgemm)("T","N",&rank,&width,&n,&one,c->projected_factor,&n,
      c->projected_residual,&n,&one,coeff,&rank FCONE FCONE);
    for(int f=0;f<width;++f) {
    if((f & 15)==0)np_progress_bandwidth_loop_step();
    long double residual=0,norm=0;
    for(int i=0;i<n;++i) {
      residual+=np_gnn_projected_residual(c->projected_factor+i,n,
        coeff+(size_t)f*rank,rank,c->right[(size_t)f*n+i]);
      norm+=fabs(c->right[(size_t)f*n+i]);
    }
    const long double bound=integral_bound*residual*(2*norm+residual)*
      (1+(4*n+16)*DBL_EPSILON);
    max_bound=fmaxl(max_bound,bound);
    q->norm_max=fmax(q->norm_max,nextafter((double)(norm+residual)*(1+32*DBL_EPSILON),R_PosInf));
    }
  }

  /* Reserve at most 3/4 of the unchanged 1e-12 total
   * absolute budget for the certified compression error. */
  if(!R_FINITE((double)max_bound))return 1;
  if(max_bound>7.5e-13L) {

    return 2;
  }
  q->absolute_scale=.25;
  q->compression_bound=(double)max_bound;
  q->rank=rank;

  return 0;
}

static void np_gnn_projected_multiply(NPGNNProjectedIntegral *q,
  int nq,int points,int offset,int begin,int width)
{
  NPGNNConditionalProjectionOwner *c=q->owner;
  const int n=c->n,output_stride=nq*q->groups;
  const double one=1,zero=0;
  if(!width)return;
  if(q->rank) {
    F77_CALL(dgemm)("N","N",&points,&width,&q->rank,&one,c->profile_right+offset,&nq,
      c->profile_left+(size_t)(q->first+begin)*q->rank,&q->rank,&zero,
      c->products+(size_t)begin*output_stride+offset,&output_stride FCONE FCONE);
  } else for(int g=0;g<q->groups;++g) {
    const int start=c->projected_start[g],donors=c->projected_start[g+1]-start;
    F77_CALL(dgemm)("N","N",&points,&width,&donors,&one,c->overlap+(size_t)start*nq+offset,&nq,
      c->right+(size_t)begin*n+start,&n,&zero,
      c->products+(size_t)begin*output_stride+g*nq+offset,&output_stride FCONE FCONE);
  }
  if(q->rank)
    for(int f=begin;f<begin+width;++f)for(int z=offset;z<offset+points;++z)
      c->products[(size_t)f*output_stride+z]-=
        c->projected_diagonal[q->first+f]*
          c->overlap[(size_t)c->geometry->order[q->first+f]*nq+z];
}

static double np_gnn_projected_square(
  const NPGNNProjectedIntegral *q,const double *v,int node,int stride)
{
  if(q->groups==1 && !q->owner->has_categories)return v[node]*v[node];
  long double value=0;
  for(int g=0;g<q->groups;++g)for(int h=0;h<=g;++h)
    value+=(g==h?1:2)*q->owner->profile_moments[(size_t)g*q->groups+h]*
      (long double)v[(size_t)g*stride+node]*v[(size_t)h*stride+node];
  return (double)value;
}

/* A fast-varying factor far from its Gaussian peak cannot
 * hide material mass between quadrature nodes. Nothing is truncated: the
 * canonical factors are still evaluated. This bounds the possible unobserved
 * integral AND quadrature contribution before waiving only the variation
 * guard; the same absolute/relative integration budget is retained. */
static double np_gnn_projected_guard(NPGNNProjectedIntegral *q,
  double lo,double hi)
{
  const NPGNNIntegralGeometry *g=q->owner->geometry;
  if(q->owner->kernel!=0) {
    /* Same signed-kernel Gauss8 remainder bound as the prefix owner.
     * L1 bounds include any certified compression residual. */
    NPGNNConditionalPrefix bound={.n=q->owner->n,.y=q->owner->data[0],
      .yscale=g->scale,.anchor=q->anchor,.maxL1=q->norm_max};
    memcpy(bound.logG,q->owner->logG,sizeof(bound.logG));
    return np_cgnn_gaussian_interval_bound(&bound,lo,hi);
  }
  const double width=hi-lo,scale=g->scale;
  const double variation=fmax(fabs(q->anchor-g->sorted[0]),
                              fabs(q->anchor-g->sorted[g->n-1]));
  if(variation*width/scale<=4)return 0;
  double log_maximum=R_NegInf;
  for(int i=0;i<g->n;++i) {
    const double slope=(q->anchor-q->owner->data[0][i])/scale;
    if(fabs(slope)*width<=4)continue;
    const double a=1/scale+slope*lo,b=1/scale+slope*hi;
    const double closest=(a<=0 && b>=0)||(b<=0 && a>=0)?0:fmin(fabs(a),fabs(b));
    log_maximum=fmax(log_maximum,log(allck[0](0))-log(scale)-.5*closest*closest);
  }
  if(q->norm_max==0 || log_maximum==R_NegInf)return 0;
  const double log_full=log(q->norm_max)+log(allck[0](0))-log(scale);
  const double log_unseen=log(q->norm_max)+log_maximum;
  const double log_bound=log(q->category_max)+log(width)+log(4.0)+log_full+log_unseen+
    log1p(.5*exp(log_unseen-log_full));
  return log_bound<log(DBL_MIN)?DBL_MIN:exp(log_bound);
}

static void np_gnn_projected_rule(NPGNNProjectedIntegral *q,
  double lo,double hi,int rule,double *out)
{
  NPGNNConditionalProjectionOwner *c=q->owner;
  const int nq=rule?8:4,n=c->n,b=q->rows;
  const double *nodes=rule?np_gnn_integral_nodes8:np_gnn_integral_nodes4;
  const double *weights=rule?np_gnn_integral_weights8:np_gnn_integral_weights4;
  const double scale=c->geometry->scale,mid=.5*lo+.5*hi,half=.5*hi-.5*lo;
  const double one=1,zero=0;
  for(int i=0;i<n;++i)for(int z=0;z<nq;++z)
    c->left[i*nq+z]=(1+(q->anchor-c->data[0][c->projected_order[i]])*(mid+half*nodes[z]))/scale;

  np_ckernelv(c->kernel,c->left,n*nq,0,0,1,c->overlap,NULL,0,0,1,1/scale,
    0,0,0,NULL,NULL);

  const int output_stride=nq*q->groups;
  memset(c->products,0,(size_t)output_stride*b*sizeof(double));
  if(q->rank)
    F77_CALL(dgemm)("N","N",&nq,&q->rank,&n,&one,c->overlap,&nq,
      c->projected_factor,&n,&zero,c->profile_right,&nq FCONE FCONE);
  for(int begin=0;begin<b;) {
    while(begin<b && !q->selected[begin])++begin;
    int end=begin;
    while(end<b && q->selected[end])++end;
    const int width=end-begin;
    np_gnn_projected_multiply(q,nq,nq,0,begin,width);
    begin=end;
  }
  for(int f=0;f<b;++f) {
    double value=0;
    if(q->selected[f])for(int z=0;z<nq;++z)
      value+=weights[z]*np_gnn_projected_square(q,
        c->products+(size_t)f*output_stride,z,nq);
    out[f]=half*value;
  }
}

static int np_gnn_projected_integrate(NPGNNProjectedIntegral *q,
  double lo,double hi,int depth,double *out)
{
  const size_t capacity=(size_t)q->owner->projected_rows;
  double *low=q->owner->projected_recursive+3*(size_t)depth*capacity;
  double *high=low+capacity,*other=high+capacity;int pass=1;
  const int certificate_domain=q->owner->kernel==0 && q->groups==1 && q->owner->filled_row != NULL &&
    q->owner->geometry->scale==1;
  const double analytic=certificate_domain ? np_cgnn_gauss4_log_bound(lo,hi,
    q->anchor-q->owner->geometry->sorted[q->owner->n-1],
    q->anchor-q->owner->geometry->sorted[0],q->norm_max*q->norm_max) : R_PosInf;
  if(certificate_domain &&
     analytic<=log(ldexp(q->absolute_scale*1e-12/(8*q->owner->n),-depth))){
    np_gnn_projected_rule(q,lo,hi,0,out);
    for(int f=0;f<q->rows;++f)if(q->selected[f]){
      if(!R_FINITE(out[f]))return 1;
      q->owner->bounds[q->owner->geometry->order[q->first+f]]+=exp(analytic);
    }
    return 0;
  }
  np_gnn_projected_rule(q,lo,hi,0,low);
  np_gnn_projected_rule(q,lo,hi,1,high);
  for(int f=0;f<q->rows;++f) {
    if(!R_FINITE(high[f]))return 1;
    if(fabs(low[f]-high[f])>ldexp(q->absolute_scale*1e-12/(4*q->owner->n),-depth)+
       1e-10*fabs(high[f]))pass=0;
  }
  const double unseen=np_gnn_projected_guard(q,lo,hi);
  const double budget=ldexp(q->absolute_scale*1e-12/(4*q->owner->n),-depth);
  const int guard_fail=!R_FINITE(unseen) || unseen>budget/4;
  for(int f=0;f<q->rows;++f)if(q->selected[f] &&
    fabs(low[f]-high[f])+unseen>budget+1e-10*fabs(high[f]))pass=0;
  if(guard_fail)pass=0;
  if(pass){
    for(int f=0;f<q->rows;++f)if(q->selected[f])
      q->owner->bounds[q->owner->geometry->order[q->first+f]]+=fabs(low[f]-high[f])+unseen;
    memcpy(out,high,(size_t)q->rows*sizeof(double));return 0;
  }
  const double mid=.5*lo+.5*hi;
  if(depth>=20 || !(lo<mid && mid<hi)) {
    int fold=0;while(fold<q->rows && !q->selected[fold])++fold;
    np_conditional_failure_record(NP_CONDITIONAL_WORK_EXHAUSTED,
      fold<q->rows ? 1+q->owner->geometry->order[q->first+fold] : 0);
    return 1;
  }
  if(np_gnn_projected_integrate(q,lo,mid,depth+1,out) ||
     np_gnn_projected_integrate(q,mid,hi,depth+1,other))return 1;
  for(int f=0;f<q->rows;++f)out[f]+=other[f];
  return 0;
}

typedef struct {
  int count;
  double origin[192],offset[192],anchor[192],weight[192];
  int first_deleted[192],end_deleted[192],successor[192];
} NPGNNProjectedCompactBatch;

/* Shared response nodes, then canonical influence contraction, then square.
 * Computing a full bounded BLAS slab also permits one batch to span geometry
 * intervals. The exact deleted-radius mask is applied only at accumulation. */
static int np_gnn_projected_compact_batch(NPGNNProjectedIntegral *q,
  NPGNNProjectedCompactBatch *b)
{
  NPGNNConditionalProjectionOwner *c=q->owner;
  const int nq=b->count,n=c->n;
  if(!nq)return 0;
  for(int i=0;i<n;++i)for(int z=0;z<nq;++z) {
    const double slope=b->anchor[z]-c->data[0][c->projected_order[i]];
    c->left[(size_t)i*nq+z]=
      fma(slope,b->offset[z],fma(slope,b->origin[z],1))/c->geometry->scale;
  }
  np_ckernelv(c->kernel,c->left,n*nq,0,0,1,c->overlap,NULL,0,0,1,
    1/c->geometry->scale,0,0,0,NULL,NULL);
  const double one=1,zero=0;
  if(q->rank)
    F77_CALL(dgemm)("N","N",&nq,&q->rank,&n,&one,c->overlap,&nq,
      c->projected_factor,&n,&zero,c->profile_right,&nq FCONE FCONE);
  np_gnn_projected_multiply(q,nq,nq,0,0,q->rows);
  for(int f=0;f<q->rows;++f) {
    const int position=q->first+f,fold=c->geometry->order[position];
    for(int z=0;z<nq;++z) {
      const int deleted=position>=b->first_deleted[z] && position<b->end_deleted[z];
      if(deleted!=b->successor[z])continue;
      const double fit=c->products[(size_t)f*nq+z];
      const double contribution=b->weight[z]*fit*fit;
      if(!R_FINITE(contribution))return 1;
      np_cgnn_compensated(contribution,c->values+fold,c->projected_recursive+f);
    }
  }
#ifdef NP_CF167_TRACE
  np_cgnn_trace_projected_nodes+=(unsigned long)nq;
#endif
  b->count=0;
  return 0;
}

/* Compact kernels have degree d on each support piece. Their squared fit
 * needs d+1 Gauss nodes, with no adaptive error decision. Support cuts are
 * shared by all folds using the same primary/successor radius. */
static int np_gnn_projected_compact(NPGNNProjectedIntegral *q)
{
  NPGNNConditionalProjectionOwner *c=q->owner;
  const NPGNNIntegralGeometry *g=c->geometry;
  const int rule=c->kernel==8 ? 1 : 2*(c->kernel-3)+1;
  const double *nodes=rule==1 ? np_cgnn_nodes1 : rule==3 ? np_cgnn_nodes3 :
    rule==5 ? np_cgnn_nodes5 : rule==7 ? np_cgnn_nodes7 : np_cgnn_nodes9;
  const double *weights=rule==1 ? np_cgnn_weights1 : rule==3 ? np_cgnn_weights3 :
    rule==5 ? np_cgnn_weights5 : rule==7 ? np_cgnn_weights7 : np_cgnn_weights9;
  const double support=c->kernel==8 ? 1 : sqrt(5.0);
  NPGNNProjectedCompactBatch batch={0};
  if(q->groups!=1 || c->has_categories || !c->folded)return 1;
  memset(c->projected_recursive,0,(size_t)q->rows*sizeof(double));
  for(int at=0;at<g->count;++at) {
    np_progress_bandwidth_loop_step();
    const NPGNNIntegralInterval t=g->intervals[at];
    const int begin=MAX(q->first,t.first_deleted);
    const int end=MIN(q->first+q->rows,t.end_deleted);
    const int deleted=MAX(0,end-begin);
    for(int successor=0;successor<2;++successor) {
      if((successor && !deleted) || (!successor && deleted==q->rows))continue;
      const double anchor=successor ? t.successor_anchor : t.primary_anchor;
      double lo,hi,width;
      if(np_gnn_integral_reciprocal_interval(t.lo,t.hi,anchor,&lo,&hi,&width) ||
         !R_FINITE(lo) || !R_FINITE(hi))return 1;
      if(width==0)continue;
      int cuts=0;
      c->compact_cuts[cuts++]=0;
      c->compact_cuts[cuts++]=width;
      for(int i=0;i<c->n;++i) {
        const double slope=anchor-c->data[0][i];
        if(slope==0)continue;
        for(int sign=-1;sign<=1;sign+=2) {
          /* Store an offset from lo. This also retains intervals whose
           * separately rounded reciprocals coincide but width is positive. */
          const double offset=fma(-slope,lo,sign*support*g->scale-1)/slope;
          if(0<offset && offset<width)c->compact_cuts[cuts++]=offset;
        }
      }
      qsort(c->compact_cuts,cuts,sizeof(double),np_cgnn_compare_cut);
      for(int piece=1;piece<cuts;++piece) {
        const double left=c->compact_cuts[piece-1],right=c->compact_cuts[piece];
        if(!(left<right))continue;
        if(batch.count+rule>192 && np_gnn_projected_compact_batch(q,&batch))return 1;
        const double mid=.5*left+.5*right,half=.5*(right-left);
        for(int z=0;z<rule;++z) {
          const int slot=batch.count++;
          batch.origin[slot]=lo;batch.offset[slot]=fma(half,nodes[z],mid);
          batch.anchor[slot]=anchor;batch.weight[slot]=half*weights[z];
          batch.first_deleted[slot]=t.first_deleted;
          batch.end_deleted[slot]=t.end_deleted;
          batch.successor[slot]=successor;
        }
      }
    }
  }
  if(np_gnn_projected_compact_batch(q,&batch))return 1;
  for(int f=0;f<q->rows;++f) {
    const int fold=g->order[q->first+f];
    c->values[fold]+=c->projected_recursive[f];
    if(!R_FINITE(c->values[fold]))return 1;
  }
#ifdef NP_CF167_TRACE
  REprintf("CF167 compact projected rank=%d first_fold=%d end_fold=%d kernel=%d compression_rank=%d\n",
    np_cgnn_trace_rank(),q->first,q->first+q->rows,c->kernel,q->rank);
#endif
  return 0;
}

static int np_gnn_projected_contract(NPGNNConditionalProjectionOwner *c)
{
  const int n=c->n,folds=c->folded?n:1;
  const int groups=c->has_categories?c->nprofiles:1;
  int block=c->projected_rows;
  const int batch=c->projected_intervals,nodes=batch*24;
  if(groups<1 || block<1 || batch<1 || c->projected_bytes>c->budget)return 1;
  if(c->first_fold==c->end_fold)return 0;
  c->left=(double *)malloc((size_t)nodes*n*sizeof(double));
  c->overlap=(double *)malloc((size_t)nodes*n*sizeof(double));
  c->right=(double *)malloc((size_t)n*MIN(512,block)*sizeof(double));
  c->products=(double *)malloc((size_t)nodes*block*groups*sizeof(double));
  c->weight_row=(double *)malloc((size_t)n*sizeof(double));
  c->projected_order=(int *)malloc((size_t)n*sizeof(int));
  c->projected_start=(int *)calloc((size_t)groups+1,sizeof(int));
  c->profile_moments=(double *)malloc((size_t)groups*groups*sizeof(double));
  c->projected_selected=(int *)calloc((size_t)block,sizeof(int));
  c->projected_recursive=(double *)malloc((size_t)block*3*21*sizeof(double));
  c->projected_estimates=(double *)malloc((size_t)block*2*16*sizeof(double));
  if(c->kernel>=4) {
    c->compact_cuts=(double *)malloc((2*(size_t)n+2)*sizeof(double));
    if(!c->compact_cuts)return 1;
  }
  if(c->projected_shared) {
    c->projected_factor=(double *)malloc((size_t)n*512*sizeof(double));
    c->projected_diagonal=(double *)calloc(n,sizeof(double));
    c->profile_right=(double *)malloc((size_t)nodes*512*sizeof(double));
    c->projected_residual=(double *)malloc((size_t)n*512*sizeof(double));
    if(!c->projected_factor || !c->projected_diagonal ||
       !c->profile_right || !c->projected_residual)return 1;
  }
  if(!c->left || !c->overlap || !c->right || !c->products || !c->weight_row ||
     !c->projected_order || !c->projected_start || !c->profile_moments ||
     !c->projected_selected || !c->projected_recursive || !c->projected_estimates)return 1;
  const NPGNNIntegralGeometry *g=c->geometry;
  int position=0;
  double category_max=1;
  for(int group=0;group<groups;++group) {
    c->projected_start[group]=position;
    for(int i=0;i<n;++i)if(!c->has_categories || c->profile_id[i]==group)
      c->projected_order[position++]=i;
    for(int other=0;other<groups;++other) {
      const double value=c->has_categories?
        c->category_overlap(c->context,c->profile_rep[group],c->profile_rep[other]):1;
      if(!R_FINITE(value))return 1;
      c->profile_moments[(size_t)group*groups+other]=value;
      category_max=fmax(category_max,fabs(value));
    }
  }
  c->projected_start[groups]=position;
  if(position!=n)return 1;
  NPGNNProjectedIntegral shared={.owner=c,.rows=folds,
    .groups=groups,.category_max=category_max,.absolute_scale=1,
    .selected=c->projected_selected};
  if(c->projected_shared) {
    const int prepared=np_gnn_projected_factor(&shared);
    if(prepared==1)return 1;
    if(prepared==2) {
      for(int f=0;f<folds;++f)if(c->values[f]!=0.0)return 1;
      free(c->projected_factor);c->projected_factor=NULL;
      free(c->projected_diagonal);c->projected_diagonal=NULL;
      free(c->projected_residual);c->projected_residual=NULL;
      free(c->profile_left);c->profile_left=NULL;
      free(c->profile_right);c->profile_right=NULL;
      c->projected_shared=0;block=MIN(512,folds);
    }
  }
  for(int first=c->first_fold;first<c->end_fold;first+=block) {
    NPGNNProjectedIntegral q={.owner=c,.rows=MIN(block,c->end_fold-first),
      .groups=groups,.category_max=category_max,.absolute_scale=1,
      .first=first,.selected=c->projected_selected};
    if(c->projected_shared) {
      q=shared;q.first=first;q.rows=MIN(block,c->end_fold-first);
    }
    else for(int f=0;f<q.rows;++f) {
      if(c->row(c->context,c->folded?g->order[first+f]:0,c->weight_row))return 1;
      double norm=0;
      for(int i=0;i<n;++i) {
        const double value=c->weight_row[c->projected_order[i]];
        if(!R_FINITE(value))return 1;
        c->right[(size_t)f*n+i]=value;norm+=fabs(value);
      }
      q.norm_max=fmax(q.norm_max,norm);
    }

    if(c->kernel>=4) {
      if(np_gnn_projected_compact(&q))return 1;
      continue;
    }
    for(int start=0;start<g->count;start+=batch) {
      np_progress_bandwidth_loop_step();
      double left[16],right[16],anchor[16],interval_width[16],node[192],node_anchor[192];
      int interval[16],alternate[16],records=0,offset[16],points[16],packed=0;
      double analytic[16];
#ifdef NP_CF167_TRACE
      for(int at=0;at<192;++at)node[at]=node_anchor[at]=R_NaN;
#endif
      for(int z=start;z<MIN(start+batch,g->count);++z) {
        if(!c->folded && z%c->pair_ranks!=c->pair_rank)continue;
        const NPGNNIntegralInterval t=g->intervals[z];
        for(int successor=0;successor<(c->folded?2:1);++successor) {
          int any=0;
          for(int f=0;f<q.rows;++f) {
            const int pos=first+f;
            const int deleted=c->folded && pos>=t.first_deleted && pos<t.end_deleted;
            any|=deleted==successor;
          }
          if(!any)continue;
          double width;
          anchor[records]=successor?t.successor_anchor:t.primary_anchor;
          if(np_gnn_integral_reciprocal_interval(t.lo,t.hi,anchor[records],
               left+records,right+records,&width) ||
             !R_FINITE(left[records]) || !R_FINITE(right[records]))return 1;
          interval_width[records]=width;
          interval[records]=z;alternate[records]=successor;
          const double mid=.5*left[records]+.5*right[records],half=.5*width;
          const int certificate_domain=c->kernel==0 && q.groups==1 && c->filled_row != NULL && g->scale==1;
          analytic[records]=certificate_domain ? np_cgnn_gauss4_log_bound(left[records],right[records],
             anchor[records]-g->sorted[n-1],anchor[records]-g->sorted[0],q.norm_max*q.norm_max) : R_PosInf;
          const int certified=certificate_domain &&
            analytic[records]<=log(q.absolute_scale*1e-12/(8*n));
          offset[records]=packed;points[records]=certified?4:12;
          for(int zq=0;zq<points[records];++zq) {
            node[packed+zq]=mid+half*(zq<4?np_gnn_integral_nodes4[zq]:np_gnn_integral_nodes8[zq-4]);
            node_anchor[packed+zq]=anchor[records];
          }
          packed+=points[records];
          ++records;
        }
      }
      const int nq=packed;
      if(nq==0)continue;
#ifdef NP_CF167_TRACE
      if(nq>192)return 1;
      for(int at=0;at<nq;++at)if(!R_FINITE(node[at]) || !R_FINITE(node_anchor[at]))
        error("diagnostic GNN node packing left an undefined slot");
      np_cgnn_trace_projected_nodes+=(unsigned long)nq;
      for(int rec=0;rec<records;++rec)np_cgnn_trace_certificates+=points[rec]==4;
#endif
      for(int i=0;i<n;++i)for(int zq=0;zq<nq;++zq)
        c->left[i*nq+zq]=(1+(node_anchor[zq]-c->data[0][c->projected_order[i]])*node[zq])/g->scale;

      np_ckernelv(c->kernel,c->left,n*nq,0,0,1,c->overlap,NULL,0,0,1,1/g->scale,
        0,0,0,NULL,NULL);

      const double one=1,zero=0;
      const int output_stride=nq*groups;
      memset(c->products,0,(size_t)output_stride*q.rows*sizeof(double));
      if(q.rank)
        F77_CALL(dgemm)("N","N",&nq,&q.rank,&n,&one,c->overlap,&nq,
          c->projected_factor,&n,&zero,c->profile_right,&nq FCONE FCONE);
      for(int rec=0;rec<records;++rec) {
        const NPGNNIntegralInterval t=g->intervals[interval[rec]];
        const int active_begin=MAX(0,MIN(q.rows,t.first_deleted-first));
        const int active_end=MAX(0,MIN(q.rows,t.end_deleted-first));
        for(int part=0;part<(alternate[rec]?1:2);++part) {
          const int begin=alternate[rec]?active_begin:(part?active_end:0);
          const int end=alternate[rec]?active_end:(part?q.rows:active_begin);
          const int width=end-begin;
          np_gnn_projected_multiply(&q,nq,points[rec],offset[rec],begin,width);
        }
      }
      /* Copy each projected slab: recursive refinement reuses products. */
      double *estimates=c->projected_estimates;
      double *errors=estimates+(size_t)16*c->projected_rows;
      for(int rec=0;rec<records;++rec)for(int f=0;f<q.rows;++f) {
        const double *v=c->products+(size_t)f*output_stride+offset[rec];
        double low=0,high=0;
        for(int zq=0;zq<4;++zq)low+=np_gnn_integral_weights4[zq]*
          np_gnn_projected_square(&q,v,zq,nq);
        if(points[rec]==4)high=low;
        else for(int zq=0;zq<8;++zq)high+=np_gnn_integral_weights8[zq]*
          np_gnn_projected_square(&q,v,4+zq,nq);
        estimates[(size_t)rec*q.rows+f]=high*interval_width[rec]*.5;
        errors[(size_t)rec*q.rows+f]=points[rec]==4?exp(analytic[rec]):fabs(high-low)*interval_width[rec]*.5;
      }
      for(int rec=0;rec<records;++rec) {
        const NPGNNIntegralInterval t=g->intervals[interval[rec]];
        int pass=1;
        q.anchor=anchor[rec];
        for(int f=0;f<q.rows;++f) {
          const int pos=first+f;
          const int deleted=c->folded && pos>=t.first_deleted && pos<t.end_deleted;
          q.selected[f]=deleted==alternate[rec];
          if(!R_FINITE(estimates[(size_t)rec*q.rows+f]))return 1;
          if(q.selected[f] && errors[(size_t)rec*q.rows+f]>q.absolute_scale*1e-12/(4*n)+
             1e-10*fabs(estimates[(size_t)rec*q.rows+f]))pass=0;
        }
        const double unseen=points[rec]==4?0:np_gnn_projected_guard(&q,left[rec],right[rec]);
        const int guard_fail=!R_FINITE(unseen) || unseen>q.absolute_scale*1e-12/(16*n);
        for(int f=0;f<q.rows;++f)if(q.selected[f] &&
          errors[(size_t)rec*q.rows+f]+unseen>q.absolute_scale*1e-12/(4*n)+
            1e-10*fabs(estimates[(size_t)rec*q.rows+f]))pass=0;
        if(guard_fail)pass=0;
        if(!pass) {
          if(np_gnn_projected_integrate(&q,left[rec],right[rec],0,
               estimates+(size_t)rec*q.rows))return 1;
        }
        if(pass)for(int f=0;f<q.rows;++f)if(q.selected[f])
          c->bounds[g->order[first+f]]+=errors[(size_t)rec*q.rows+f]+unseen;
        for(int f=0;f<q.rows;++f)if(q.selected[f])
          c->values[c->folded?g->order[first+f]:0]+=estimates[(size_t)rec*q.rows+f];
      }
    }
  }
  if(c->projected_shared)
    for(int f=c->first_fold;f<c->end_fold;++f)
      c->bounds[g->order[f]]+=shared.compression_bound;

  if(c->has_categories && !c->folded)c->pair_rows[0]=c->values[0];
  return 0;
}


typedef struct {
  NPConditionalCVLSRouteContext *route;
  NPGNNIntegralUniformContext categories;
  double *values,*bounds,*xrow,*yrow,*cv;
  int *support_reps;
  int status;
} NPGNNConditionalIntegralContext;

/* Structural admission for this conditional-density GNN CVLS owner only.
 * The incumbent provider leaves its actual X weights and basis in X-tree
 * order. Count distinct complete deleted basis rows, stopping at basis width.
 * This is a rank upper bound, never a numerical full-rank certificate. */
static int np_cgnn_integral_x_row(NPGNNConditionalIntegralContext *c,
                                 int evaluation)
{
  if(np_conditional_cvls_provider_x_row(c->route,evaluation,c->xrow))return 1;
  if(BANDWIDTH_den_extern!=BW_GEN_NN || c->route->beta_x ||
     num_reg_continuous_extern<2 || int_cxker_bound_extern ||
     (KERNEL_reg_extern!=4 && KERNEL_reg_extern!=8) ||
     np_lp_engine_extern!=NP_LP_ENGINE_GENERAL)return 0;
  NPConditionalXRowCtx *x=&c->route->legacy_x;
  const int p=np_glp_cv_cache.nterms, n=num_obs_train_extern;
  const int held=int_TREE_X==NP_TREE_TRUE ? ipt_lookup_extern_X[evaluation] : evaluation;
  int count=0;
  if(p<=1)return 0;
  if(!c->support_reps)c->support_reps=np_cgnn_calloc(p,sizeof(int));
  if(!c->support_reps)return 1;
  for(int j=0;j<n && count<p;++j) {
    if(j==held || x->kw[j]==0.0)continue;
    int duplicate=0;
    for(int a=0;a<count && !duplicate;++a) {
      const int other=c->support_reps[a];
      int k=0;
      while(k<p && x->basis[k][j]==x->basis[k][other])++k;
      duplicate=k==p;
    }
    if(!duplicate)c->support_reps[count++]=j;
  }
  return count<p;
}

/* The route produces original-index X rows; the integral has Y-tree donor
 * order. The fold and both donor maps are explicit and applied exactly once. */
static int np_gnn_conditional_integral_weight(void *raw,int fold,double *row)
{
  NPGNNConditionalIntegralContext *c=(NPGNNConditionalIntegralContext *)raw;
  const int i=int_TREE_Y == NP_TREE_TRUE ? ipt_extern_Y[fold] : fold;
  if(np_cgnn_integral_x_row(c,i))return 1;
  for(int j=0;j<num_obs_train_extern;++j)
    row[j]=c->xrow[int_TREE_Y == NP_TREE_TRUE ? ipt_extern_Y[j] : j];
  return 0;
}

/* A compression-only receipt, without widening the shared X-row interface.
 * The incumbent deleted row and normalizing order are left untouched. */
typedef struct {
  NPConditionalXRowCtx *x;
  NPConditionalBoundState bounds;
  int status;
  double self;
} NPGNNFilledSelf;

static void np_gnn_filled_self_cleanup(void *raw,Rboolean jump)
{
  NPGNNFilledSelf *c=(NPGNNFilledSelf *)raw;
  (void)jump;
  np_conditional_pop_bounds(&c->bounds);
}

static SEXP np_gnn_filled_self_execute(void *raw)
{
  NPGNNFilledSelf *c=(NPGNNFilledSelf *)raw;
  NPConditionalXRowCtx *x=c->x;
  c->status=np_conditional_kernel_row_raw(
    x->kernel_cx,x->kernel_ux,x->kernel_ox,x->x_operator,BW_GEN_NN,1,
    num_reg_unordered_extern,num_reg_ordered_extern,num_reg_continuous_extern,
    x->eval_xuno_one,x->eval_xord_one,x->eval_xcon_one,
    x->eval_xuno_one,x->eval_xord_one,x->eval_xcon_one,
    x->vsfx,1,x->matrix_bandwidth_eval_one,x->matrix_bandwidth_eval_one,
    x->lambdax,num_categories_extern_X,matrix_categorical_vals_extern_X,
    NP_TREE_FALSE,NULL,&c->self,NULL);
  return R_NilValue;
}

static int np_gnn_conditional_integral_filled_weight(void *raw,int fold,
  double *row,double *diagonal)
{
  NPGNNConditionalIntegralContext *c=(NPGNNConditionalIntegralContext *)raw;
  if(c->route==NULL || !c->route->ready || c->route->beta_x ||
     BANDWIDTH_den_extern!=BW_GEN_NN || np_lp_engine_extern!=NP_LP_ENGINE_SCALAR)
    return 1;
  if(np_gnn_conditional_integral_weight(raw,fold,row))return 1;
  NPConditionalXRowCtx *x=&c->route->legacy_x;
  if(x->num_reg_tot<=0) {
    *diagonal=1.0/(num_obs_train_extern-1);
    return 0;
  }
  double sum=0;
  for(int i=0;i<num_obs_train_extern;++i)sum+=x->kw[i];
  if(!(fabs(sum)>DBL_MIN))return 1;
  NPGNNFilledSelf self={.x=x,.status=1,.self=R_NaN};
  np_conditional_push_bounds(int_cxker_bound_extern,
    vector_cxkerlb_extern,vector_cxkerub_extern,&self.bounds);
  R_UnwindProtect(np_gnn_filled_self_execute,&self,
    np_gnn_filled_self_cleanup,&self,NULL);
  if(self.status || !R_FINITE(self.self))return 1;
  *diagonal=self.self/sum;
  return !R_FINITE(*diagonal);
}

typedef struct {
  NPGNNConditionalIntegralContext adapter;
  NPGNNConditionalProjectionOwner integral;
  NPGNNIntegralGeometry geometry;
  double *cross,score;
  int status;
} NPGNNConditionalProjectedCall;

static void np_cgnn_projected_cleanup(void *raw,Rboolean jump)
{
  NPGNNConditionalProjectedCall *a=raw;
  NPGNNConditionalProjectionOwner *c=&a->integral;
  (void)jump;
  np_gnn_integral_geometry_clear(&a->geometry);
  free(a->adapter.xrow);free(a->adapter.yrow);free(a->cross);
  free(a->adapter.support_reps);
  free(c->values);free(c->bounds);free(c->weight_row);
  free(c->left);free(c->right);free(c->overlap);free(c->products);
  free(c->profile_moments);free(c->profile_left);free(c->profile_right);
  free(c->projected_order);free(c->projected_start);
  free(c->projected_selected);free(c->projected_recursive);
  free(c->projected_estimates);free(c->projected_factor);
  free(c->projected_diagonal);free(c->projected_residual);
  free(c->compact_cuts);
}

static int np_cgnn_projected_local(NPGNNConditionalProjectedCall *a)
{
  NPGNNConditionalProjectionOwner *c=&a->integral;
  NPGNNConditionalIntegralContext *b=&a->adapter;
  const int n=num_obs_train_extern;
  *c=(NPGNNConditionalProjectionOwner){.n=n,.folded=1,.first_fold=0,.kernel=KERNEL_den_extern,
    .end_fold=n,.pair_ranks=1,.data=matrix_Y_continuous_train_extern,
    .row=np_gnn_conditional_integral_weight,
    .filled_row=!b->route->beta_x && np_lp_engine_extern==NP_LP_ENGINE_SCALAR ?
      np_gnn_conditional_integral_filled_weight:NULL,
    .geometry=&a->geometry,.context=b,.budget=NP_CONDITIONAL_LP_TILE_BUDGET_BYTES};
  if(c->kernel>0 && c->kernel<4)np_cgnn_gaussian_derivative_constants(c->kernel,c->logG);
  c->values=np_cgnn_calloc(n,sizeof(double));
  c->bounds=np_cgnn_calloc(n,sizeof(double));
  b->xrow=np_cgnn_calloc(n,sizeof(double));
  b->yrow=np_cgnn_calloc(n,sizeof(double));
  a->cross=np_cgnn_calloc(n,sizeof(double));
  if(!c->values || !c->bounds || !b->xrow || !b->yrow || !a->cross)return 1;
  if(np_gnn_integral_geometry_prepare(&a->geometry,c->data[0],n,
       b->route->legacy_y.vsfy[0],1))return 1;
  /* Finite reciprocal panels admit the projected representation. Compact
   * kernels can have finite integrals even with zero-radius endpoints
   * (notably uniform k=1/ties); retain their existing whole-support pair
   * owner before any projected row or integral is computed. */
  for(int z=0;z<a->geometry.count;++z) {
    const NPGNNIntegralInterval t=a->geometry.intervals[z];
    if(t.lo==t.primary_anchor || t.hi==t.primary_anchor ||
       t.lo==t.successor_anchor || t.hi==t.successor_anchor) {
      if(c->kernel>=4) {
        c->compact_pair_geometry=1;
        return 0;
      }
      return 1;
    }
  }
  if(!np_gnn_projected_storage(c))return 1;
#ifdef MPI2
  int count;
  np_objective_outer_owned_rows(0,n,np_objective_outer_rows_enabled(1),
    &c->first_fold,&count);
  c->end_fold=c->first_fold+count;
#endif
#ifdef NP_CF167_TRACE
  REprintf("CF167 projected rank=%d first_fold=%d end_fold=%d n=%d\n",
    np_cgnn_trace_rank(),c->first_fold,c->end_fold,n);
#endif
  return np_gnn_projected_contract(c);
}

static int np_cgnn_general_cvls(NPConditionalCVLSRouteContext *route,double *cv);

static SEXP np_cgnn_projected_body(void *raw)
{
  NPGNNConditionalProjectedCall *a=raw;
  NPGNNConditionalProjectionOwner *c=&a->integral;
  NPGNNConditionalIntegralContext *b=&a->adapter;
  const int n=num_obs_train_extern;
  int fail=np_cgnn_projected_local(a),first=0,count=n;
#ifdef MPI2
  const int parallel=np_objective_outer_rows_enabled(1);
  if(np_conditional_outer_buffer_finish(parallel,n,fail,c->values,NULL,
       "conditional GNN projected integral"))return R_NilValue;
  if(parallel) {
    np_mpi_allreduce_in_place_double(c->bounds,n,MPI_SUM,
      "conditional GNN projected error estimate");
    np_mpi_allreduce_in_place_double(&c->capped,1,MPI_SUM,
      "conditional GNN projected finite caps");
  }
  np_objective_outer_owned_rows(0,n,parallel,&first,&count);
#else
  if(fail)return R_NilValue;
#endif
  if(c->compact_pair_geometry) {
    a->status=np_cgnn_general_cvls(b->route,&a->score);
    return R_NilValue;
  }
  /* I2 remains the canonical deleted X/Y row dot product. */
  for(int i=first;i<first+count && !fail;++i) {
    np_progress_bandwidth_loop_step();
    double logscale=0.0,linear;
    if(np_cgnn_integral_x_row(b,i) ||
       np_conditional_cvls_provider_y_train_row(b->route,i,b->yrow,&logscale)) {
      fail=1;break;
    }
    linear=np_blas_ddot_int(n,b->xrow,b->yrow);
    if(np_continuous_kernel_scaled_restore(linear,logscale,1,&linear)!=
       NP_CONTINUOUS_ROW_OK) {fail=1;break;}
    a->cross[i]=linear;
  }
#ifdef MPI2
  if(np_conditional_outer_buffer_finish(parallel,n,fail,a->cross,NULL,
       "conditional GNN projected canonical I2"))return R_NilValue;
#else
  if(fail)return R_NilValue;
#endif
  double total=0,correction=0,estimate=0;
  for(int i=0;i<n;++i) {
    const int yi=int_TREE_Y==NP_TREE_TRUE ? ipt_lookup_extern_Y[i]:i;
    np_cgnn_compensated((c->values[yi]-2*a->cross[i])/n,&total,&correction);
    estimate+=c->bounds[yi]/n;
  }
  a->score=total+correction;
  if(!R_FINITE(a->score) || !R_FINITE(estimate))return R_NilValue;
  np_cgnn_notice_caps+=c->capped;
  if(c->capped>0)np_cgnn_notice_error=fmax(np_cgnn_notice_error,estimate);
  a->status=0;
  return R_NilValue;
}

static int np_cgnn_projected_cvls(NPConditionalCVLSRouteContext *route,double *cv)
{
  NPGNNConditionalProjectedCall call={.adapter={.route=route},.status=1};
  R_UnwindProtect(np_cgnn_projected_body,&call,np_cgnn_projected_cleanup,&call,NULL);
  if(!call.status)*cv=call.score;
  return call.status;
}
