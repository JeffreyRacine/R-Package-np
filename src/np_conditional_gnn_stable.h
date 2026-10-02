/* Conditional positive-kernel GNN local polynomial moments. The separate
 * sibling retains the existing signed-kernel owner and its storage shape. */
#include "conditional_gnn_coefficients.h"
/* Raw-coordinate local moments for the conditional GNN CVLS X contraction.
 * Included privately by np_conditional_gnn_prefix.h. This does not construct
 * a basis, solve a system, choose bandwidths, or select an estimator. */

typedef struct {
  int node;
} NPGNNStableQuery;
typedef struct {
  int n, degree, base, nodes, width;
  NPConditionalQRAcc *coefficients;
  size_t used, capacity;
  double *center, *scale, *translation, *moment, *fit, *error;
  int *ids, *ranks, *slots;
  size_t *start, *stamp, epoch;
  NPGNNStableQuery *query;

} NPGNNStableMoments;

static void *np_cgnn_calloc(size_t n, size_t width);
static double np_cgnn_stable_choose(int n, int k) {
  double value = 1.0;
  for (int j = 1; j <= k; ++j) value *= (double)(n + 1 - j) / j;
  return value;
}
static double np_cgnn_stable_power(double x, int n) {
  double value = 1.0;
  for (int j = 0; j < n; ++j) value *= x;
  return value;
}
static void np_cgnn_stable_clear(NPGNNStableMoments *c) {
  free(c->coefficients);free(c->center); free(c->scale); free(c->translation); free(c->moment);
  free(c->fit); free(c->error); free(c->ids); free(c->ranks); free(c->slots);
  free(c->start); free(c->stamp); free(c->query);
  memset(c, 0, sizeof(*c));
}
static int np_cgnn_stable_append(NPGNNStableMoments *c, int node,
                                 double x, double h) {
  if (c->used == c->capacity) return 1;
  NPGNNStableQuery *q = &c->query[c->used++];
  q->node = node;
  (void)x;(void)h;
  return 0;
}
static int np_cgnn_stable_range(NPGNNStableMoments *c, int lo, int hi,
                                double x, double h) {
  for (lo += c->base, hi += c->base; lo < hi; lo /= 2, hi /= 2) {
    if ((lo & 1) && np_cgnn_stable_append(c, lo++, x, h)) return 1;
    if ((hi & 1) && np_cgnn_stable_append(c, --hi, x, h)) return 1;
  }
  return 0;
}
static int np_cgnn_stable_prepare(NPGNNStableMoments *c, int n, int degree,
    const NPGNNConditionalOrder *order, const double *radius,
    const int *first, const int *end,const double *local_radius) {
  const int width = degree+1;
  if(degree<1 || degree>INT_MAX-1)return 1;
  c->width=width;
  int depth = 0;
  size_t count;
  c->n = n; c->degree = degree; c->base = 1;
  while (c->base < n) {
    if (c->base > INT_MAX / 4) return 1;
    c->base *= 2; ++depth;
  }
  if (c->base > INT_MAX / 2) return 1;
  c->nodes = 2 * c->base;
#define NP_CGNN_LOCAL_ALLOC(ptr, number) do { \
    c->ptr = np_cgnn_calloc((number), sizeof(*c->ptr)); \
    if (!c->ptr) return 1; \
  } while (0)
  NP_CGNN_LOCAL_ALLOC(center, c->nodes);
  NP_CGNN_LOCAL_ALLOC(scale, c->nodes);
  if (!np_size_mul_checked(c->nodes, (size_t)width * width, &count)) return 1;
  NP_CGNN_LOCAL_ALLOC(translation, count);
  if (!np_size_mul_checked(c->nodes, (size_t)width * 9, &count)) return 1;
  NP_CGNN_LOCAL_ALLOC(moment, count);
  if (!np_size_mul_checked(n, 9, &count)) return 1;
  NP_CGNN_LOCAL_ALLOC(fit, count); NP_CGNN_LOCAL_ALLOC(error, count);
  NP_CGNN_LOCAL_ALLOC(ids, n); NP_CGNN_LOCAL_ALLOC(ranks, n);
  NP_CGNN_LOCAL_ALLOC(slots, n); NP_CGNN_LOCAL_ALLOC(start, (size_t)n + 1);
  NP_CGNN_LOCAL_ALLOC(stamp, c->nodes);
  if (!np_size_mul_checked(n, (size_t)4 * (depth + 1), &c->capacity)) return 1;
  NP_CGNN_LOCAL_ALLOC(query, c->capacity);
#undef NP_CGNN_LOCAL_ALLOC
  /* Build exact raw endpoint envelopes. Padded leaves have zero weight. */
  for (int i = 0; i < c->base; ++i) {
    c->center[c->base+i] = c->scale[c->base+i] = order[i < n ? i : n-1].x;
    if (i < n) c->ids[i] = order[i].id;
  }
  for (int node = c->base-1; node; --node) {
    c->center[node] = c->center[2*node];
    c->scale[node] = c->scale[2*node+1];
  }
  for (int node = 1; node < c->nodes; ++node) {
    const double lo = c->center[node], hi = c->scale[node], mid = .5*lo + .5*hi;
    c->center[node] = mid;
    c->scale[node] = fmax(fabs(lo-mid), fabs(hi-mid));
  }
  for (int node = 2; node < c->nodes; ++node) {
    const double s = c->scale[node/2],
      a = s ? (c->center[node]-c->center[node/2])/s : 0.0,
      b = s ? c->scale[node]/s : 0.0;
    for (int p = 0; p <= degree; ++p) for (int j = 0; j <= p; ++j)
      c->translation[(size_t)node*width*width + p*width+j] =
        np_cgnn_stable_choose(p,j)*np_cgnn_stable_power(a,p-j)*np_cgnn_stable_power(b,j);
  }
  /* Recipes are indexed in sorted rank, not original identity. */
  for (int r = 0; r < n; ++r) {
    const int i = order[r].id;
    c->start[r] = c->used;
    if(!(local_radius[i]>0))continue;
    const double h = radius[i];
    if (!(h > 0.0) || !R_FINITE(h) || first[i] > r || end[i] <= r) return 1;
    c->start[r] = c->used;
    if (np_cgnn_stable_range(c, first[i], r, order[r].x, h) ||
        np_cgnn_stable_range(c, r+1, end[i], order[r].x, h)) return 1;
  }
  c->start[n] = c->used;
  c->coefficients=np_cgnn_calloc(c->used,(size_t)c->width*sizeof(*c->coefficients));
  if(!c->coefficients)return 1;
  return 0;
}

/* Visit only ancestors of the retained nonzero response donors. Empty
 * subtrees do no moment arithmetic; epoch tags exclude their stale storage. */
static void np_cgnn_stable_build(NPGNNStableMoments *c, int node, int lo, int hi,
    int first, int end, const double *basis, const double *factor, int nq) {
  if (first == end) return;
  c->stamp[node] = c->epoch;
  double *out = c->moment + (size_t)node*c->width*9;
  if (hi-lo == 1) {
    const int id = c->ids[lo], slot = c->slots[first];
    for (int z = 0; z < nq; ++z) out[z] = basis[id]*factor[slot*nq+z];
    for (int p = 1; p <= c->degree; ++p)
      for (int z = 0; z < nq; ++z) out[p*nq+z] = 0.0;
    return;
  }
  const int mid = lo+(hi-lo)/2;
  int left = first, right = end;
  while (left < right) {
    const int at = left+(right-left)/2;
    if (c->ranks[at] < mid) left = at+1; else right = at;
  }
  np_cgnn_stable_build(c, 2*node, lo, mid, first, left, basis, factor, nq);
  np_cgnn_stable_build(c, 2*node+1, mid, hi, left, end, basis, factor, nq);
  for (int p = 0; p <= c->degree; ++p) for (int z = 0; z < nq; ++z) {
    double sum = 0.0, error = 0.0;
    for (int child = 2*node; child <= 2*node+1; ++child) {
      if (c->stamp[child] != c->epoch) continue;
      const double *t = c->translation + (size_t)child*c->width*c->width;
      const double *m = c->moment + (size_t)child*c->width*9;
      for (int j = 0; j <= p; ++j)
        np_cgnn_compensated(t[p*c->width+j]*m[j*nq+z], &sum, &error);
    }
    out[p*nq+z] = sum+error;
  }
}

static NPConditionalQRAcc npc_gnn_power(NPConditionalQRAcc a,int n){
 NPConditionalQRAcc v={1,0};while(n--)v=cqr_dd_mul(v,a);return v;
}
static int npc_gnn_recipes(NPGNNStableMoments *c,int n,int k,
 const NPGNNConditionalOrder *order,const double *radius,const double *h,
 const NPConditionalQRAcc *coef,int epan,double k0){
 NPConditionalQRAcc *poly=np_cgnn_calloc(c->width,sizeof(*poly));if(!poly)return 1;
 for(int r=0;r<n;r++){
  const int i=order[r].id;if(!(radius[i]>0))continue;
  memset(poly,0,c->width*sizeof(*poly));
  for(int t=0;t<k;t++){
   NPConditionalQRAcc v=cqr_dd_mul(coef[(size_t)i*k+t],(NPConditionalQRAcc){k0,0});
   np_cqr_add(poly+t,v.hi);np_cqr_add(poly+t,v.lo);
   if(epan){double ratio=radius[i]/h[i];v=cqr_dd_mul(v,(NPConditionalQRAcc){-.2*ratio*ratio,0});np_cqr_add(poly+t+2,v.hi);np_cqr_add(poly+t+2,v.lo);}
  }
  for(size_t at=c->start[r];at<c->start[r+1];at++){
   NPGNNStableQuery *q=c->query+at;
   NPConditionalQRAcc a={c->center[q->node],0};np_cqr_add(&a,-order[r].x);a=cqr_dd_div(a,radius[i]);
   NPConditionalQRAcc b=cqr_dd_div((NPConditionalQRAcc){c->scale[q->node],0},radius[i]);
   for(int j=0;j<=c->degree;j++){
    NPConditionalQRAcc sum={0,0};
    for(int p=j;p<=c->degree;p++){
     NPConditionalQRAcc v=cqr_dd_mul(poly[p],cqr_dd_mul(npc_gnn_power(a,p-j),npc_gnn_power(b,j)));
     v=cqr_dd_mul(v,(NPConditionalQRAcc){np_cgnn_stable_choose(p,j),0});np_cqr_add(&sum,v.hi);np_cqr_add(&sum,v.lo);
    }
    if(!R_FINITE(sum.hi)||!R_FINITE(sum.lo)){free(poly);return 1;}
    c->coefficients[at*c->width+j]=sum;
   }
  }
 }
 free(poly);return 0;
}
static void np_cgnn_stable_evaluate(NPGNNStableMoments *c,const double *ones,
 const double *factor,int nq,int retained,const int *active,const double *radius){
 if(++c->epoch==0){memset(c->stamp,0,c->nodes*sizeof(size_t));c->epoch=1;}
 np_cgnn_stable_build(c,1,0,c->base,0,retained,ones,factor,nq);
 for(int r=0;r<c->n;r++){
  const int i=c->ids[r];if(!active[i] || !(radius[i]>0))continue;
  for(int z=0;z<nq;z++){
   NPConditionalQRAcc sum={0,0};
   for(size_t at=c->start[r];at<c->start[r+1];at++){
    NPGNNStableQuery *q=c->query+at;if(c->stamp[q->node]!=c->epoch)continue;
    const double *moment=c->moment+(size_t)q->node*c->width*9;
    for(int p=0;p<=c->degree;p++){
     const NPConditionalQRAcc v=c->coefficients[at*c->width+p];
     np_cqr_product(&sum,v.hi,moment[p*nq+z]);np_cqr_product(&sum,v.lo,moment[p*nq+z]);
    }
   }
   c->fit[(size_t)i*nq+z]=sum.hi+sum.lo;
  }
 }
}
