/* Raw-coordinate local moments for the conditional GNN CVLS X contraction.
 * Included privately by np_conditional_gnn_prefix.h. This does not construct
 * a basis, solve a system, choose bandwidths, or select an estimator. */
enum { NP_CGNN_LOCAL_WIDTH = 9 };
typedef struct {
  int node;
  double coefficient[NP_CGNN_LOCAL_WIDTH];
} NPGNNLocalQuery;
typedef struct {
  int n, degree, base, nodes;
  size_t used, capacity;
  double *center, *scale, *translation, *moment, *fit, *error;
  int *ids, *ranks, *slots;
  size_t *start, *stamp, epoch;
  NPGNNLocalQuery *query;
  double kernel[NP_CGNN_LOCAL_WIDTH];
} NPGNNLocalMoments;

static void *np_cgnn_calloc(size_t n, size_t width);
static double np_cgnn_local_choose(int n, int k) {
  double value = 1.0;
  for (int j = 1; j <= k; ++j) value *= (double)(n + 1 - j) / j;
  return value;
}
static double np_cgnn_local_power(double x, int n) {
  double value = 1.0;
  for (int j = 0; j < n; ++j) value *= x;
  return value;
}
static void np_cgnn_local_clear(NPGNNLocalMoments *c) {
  free(c->center); free(c->scale); free(c->translation); free(c->moment);
  free(c->fit); free(c->error); free(c->ids); free(c->ranks); free(c->slots);
  free(c->start); free(c->stamp); free(c->query);
  memset(c, 0, sizeof(*c));
}
static int np_cgnn_local_append(NPGNNLocalMoments *c, int node,
                                 double x, double h) {
  if (c->used == c->capacity) return 1;
  NPGNNLocalQuery *q = &c->query[c->used++];
  q->node = node;
  const double a = (c->center[node] - x) / h, b = c->scale[node] / h;
  for (int j = 0; j <= c->degree; ++j) {
    double sum = 0.0, error = 0.0;
    for (int p = j; p <= c->degree; ++p)
      np_cgnn_compensated(c->kernel[p] * np_cgnn_local_choose(p, j) *
          np_cgnn_local_power(a, p-j) * np_cgnn_local_power(b, j), &sum, &error);
    q->coefficient[j] = sum + error;
    if (!R_FINITE(q->coefficient[j])) return 1;
  }
  return 0;
}
static int np_cgnn_local_range(NPGNNLocalMoments *c, int lo, int hi,
                                double x, double h) {
  for (lo += c->base, hi += c->base; lo < hi; lo /= 2, hi /= 2) {
    if ((lo & 1) && np_cgnn_local_append(c, lo++, x, h)) return 1;
    if ((hi & 1) && np_cgnn_local_append(c, --hi, x, h)) return 1;
  }
  return 0;
}
static int np_cgnn_local_prepare(NPGNNLocalMoments *c, int n, int degree,
    const NPGNNConditionalOrder *order, const double *radius,
    const int *first, const int *end) {
  const int width = NP_CGNN_LOCAL_WIDTH;
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
  if (!np_size_mul_checked(c->nodes, width * width, &count)) return 1;
  NP_CGNN_LOCAL_ALLOC(translation, count);
  if (!np_size_mul_checked(c->nodes, width * 9, &count)) return 1;
  NP_CGNN_LOCAL_ALLOC(moment, count);
  if (!np_size_mul_checked(n, 9, &count)) return 1;
  NP_CGNN_LOCAL_ALLOC(fit, count); NP_CGNN_LOCAL_ALLOC(error, count);
  NP_CGNN_LOCAL_ALLOC(ids, n); NP_CGNN_LOCAL_ALLOC(ranks, n);
  NP_CGNN_LOCAL_ALLOC(slots, n); NP_CGNN_LOCAL_ALLOC(start, (size_t)n + 1);
  NP_CGNN_LOCAL_ALLOC(stamp, c->nodes);
  if (!np_size_mul_checked(n, (size_t)4 * (depth + 1), &c->capacity)) return 1;
  NP_CGNN_LOCAL_ALLOC(query, c->capacity);
#undef NP_CGNN_LOCAL_ALLOC
  if(degree==0)c->kernel[0]=.5;
  else if (degree == 2) {
    c->kernel[0] = .33541019662496845446;
    c->kernel[2] = -.067082039324993690892;
  } else {
    const double core4[] = {-15.0, 7.0}, core6[] = {2.734375, -3.28125, .721875},
      core8[] = {3.5888671875, -7.8955078125, 4.1056640625, -.5865234375};
    const double *core = degree == 4 ? core4 : degree == 6 ? core6 : core8;
    const double k = degree == 4 ? .008385254916 : .33541019662496845446;
    for (int p = 0; p < degree / 2; ++p) {
      c->kernel[2*p] += k * core[p] * (degree == 4 ? -5.0 : 1.0);
      c->kernel[2*p+2] += k * core[p] * (degree == 4 ? 1.0 : -.2);
    }
  }
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
        np_cgnn_local_choose(p,j)*np_cgnn_local_power(a,p-j)*np_cgnn_local_power(b,j);
  }
  /* Recipes are indexed in sorted rank, not original identity. */
  for (int r = 0; r < n; ++r) {
    const int i = order[r].id;
    const double h = radius[i];
    if (!(h > 0.0) || !R_FINITE(h) || first[i] > r || end[i] <= r) return 1;
    c->start[r] = c->used;
    if (np_cgnn_local_range(c, first[i], r, order[r].x, h) ||
        np_cgnn_local_range(c, r+1, end[i], order[r].x, h)) return 1;
  }
  c->start[n] = c->used;
  return 0;
}

/* Visit only ancestors of the retained nonzero response donors. Empty
 * subtrees do no moment arithmetic; epoch tags exclude their stale storage. */
static void np_cgnn_local_build(NPGNNLocalMoments *c, int node, int lo, int hi,
    int first, int end, const double *basis, const double *factor, int nq) {
  if (first == end) return;
  c->stamp[node] = c->epoch;
  double *out = c->moment + (size_t)node*NP_CGNN_LOCAL_WIDTH*9;
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
  np_cgnn_local_build(c, 2*node, lo, mid, first, left, basis, factor, nq);
  np_cgnn_local_build(c, 2*node+1, mid, hi, left, end, basis, factor, nq);
  for (int p = 0; p <= c->degree; ++p) for (int z = 0; z < nq; ++z) {
    double sum = 0.0, error = 0.0;
    for (int child = 2*node; child <= 2*node+1; ++child) {
      if (c->stamp[child] != c->epoch) continue;
      const double *t = c->translation + (size_t)child*NP_CGNN_LOCAL_WIDTH*NP_CGNN_LOCAL_WIDTH;
      const double *m = c->moment + (size_t)child*NP_CGNN_LOCAL_WIDTH*9;
      for (int j = 0; j <= p; ++j)
        np_cgnn_compensated(t[p*NP_CGNN_LOCAL_WIDTH+j]*m[j*nq+z], &sum, &error);
    }
    out[p*nq+z] = sum+error;
  }
}
static void np_cgnn_local_evaluate(NPGNNLocalMoments *c, int terms,
    const double *basis, const double *coefficient, const double *factor,
    int nq, int retained, const int *active) {
  memset(c->fit, 0, (size_t)c->n*nq*sizeof(double));
  memset(c->error, 0, (size_t)c->n*nq*sizeof(double));
  for (int t = 0; t < terms; ++t) {
    if (++c->epoch == 0) {
      memset(c->stamp, 0, (size_t)c->nodes*sizeof(size_t));
      c->epoch = 1;
    }
    np_cgnn_local_build(c, 1, 0, c->base, 0, retained, basis+(size_t)t*c->n, factor, nq);
    for (int r = 0; r < c->n; ++r) {
      const int i = c->ids[r];
      if (!active[i]) continue;
      for (int z = 0; z < nq; ++z) {
        double sum = 0.0, error = 0.0;
        for (size_t k = c->start[r]; k < c->start[r+1]; ++k) {
          const NPGNNLocalQuery *q = c->query+k;
          if (c->stamp[q->node] != c->epoch) continue;
          const double *moment = c->moment+(size_t)q->node*NP_CGNN_LOCAL_WIDTH*9;
          for (int p = 0; p <= c->degree; ++p)
            np_cgnn_compensated(q->coefficient[p]*moment[p*nq+z], &sum, &error);
        }
        const size_t at = (size_t)i*nq+z;
        np_cgnn_compensated(coefficient[(size_t)i*terms+t]*(sum+error), c->fit+at, c->error+at);
      }
    }
  }
  for (size_t j = 0; j < (size_t)c->n*nq; ++j) c->fit[j] += c->error[j];
}
