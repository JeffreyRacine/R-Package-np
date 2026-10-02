#include "np_conditional_gnn_gauss4.h"
/* Whole-support conditional GNN I1. Included only by jksum.c after the
 * canonical X/Y contexts. Local moments retain the conditioned LP basis;
 * uniform X keeps constant prefixes. Neither owner introduces a solver,
 * LP-basis conversion or response-tail cut. */
typedef struct {
  double x;
  int id;
} NPGNNConditionalOrder;
static int np_cgnn_compare(const void *a, const void *b) {
  const NPGNNConditionalOrder *x = a, *y = b;
  return (x->x > y->x) - (x->x < y->x);
}
/* Exact deleted-design upper bound for this one-X polynomial owner. Reuse
 * actual weights: signed nonzero donors count, self is excluded by identity,
 * and repeated X values cannot supply independent basis rows. Saturation is
 * not a numerical full-rank certificate. No shared solve policy is changed. */
static int np_cgnn_deleted_support_sufficient(
    const NPGNNConditionalOrder *order, const double *weights,
    const int *native_position, int n, int held, int terms) {
  int count = 0;
  double previous = 0.0;
  for (int r = 0; r < n; ++r) {
    const int id = order[r].id;
    if (id == held || weights[native_position ? native_position[id] : id] == 0.0)
      continue;
    if (count == 0 || order[r].x != previous) {
      previous = order[r].x;
      if (++count == terms)
        return 1;
    }
  }
  return 0;
}
static void np_cgnn_compensated(double v, double *s, double *e) {
  const double t = *s + v;
  *e += fabs(*s) >= fabs(v) ? (*s - t) + v : (v - t) + *s;
  *s = t;
}
#include "np_conditional_gnn_local.h"
#include "np_conditional_gnn_stable.h"
typedef struct {
  int n, folds, moments, failed, ykernel, *active, *first, *end, *orderRank;
  const double *y;
  double *x, *coef, *self, *powers;
  double maxL1, logG[17];
  NPGNNConditionalOrder *order;
  double *args, *factors, *prefix, *correction, *recursive;
  double yscale, anchor, minY, maxY, acceptedError, capped;
  int tree, selectedCount, activeCount, *coverage, *selected, *factorIndex, *activeIDs;
  int *responseRank, *responseActive, *responseSlots;
  double *cuts;
  int narrow;
  double qleft, du;
  NPGNNLocalMoments local;
  NPGNNStableMoments stable_moments;
  int stable,full_rows,deficient_rows,coefficient_count;
  size_t visits,interval_visits,work_limit;
  double *radius,*ones;NPConditionalQRAcc *stable_coef,*coefficient_work;
} NPGNNConditionalPrefix;
static void np_cgnn_local_rule(NPGNNConditionalPrefix *c, double lo, double hi,
                                int high, int compact, double *out);

/* A collapsed reciprocal panel is parameterized by t in [0,1]. Retain
 * its width separately instead of perturbing endpoints or dropping it.
 * The original q-space ratio avoids cancellation in 1+(anchor-y)*u. */
static void np_cgnn_narrow_args(NPGNNConditionalPrefix *c, int m, int nq,
                                int ordered, int pruned, double mid, double half,
                                const double *nodes) {
  for (int r = 0; r < m; ++r) {
    const int j = ordered ? c->order[pruned ? c->selected[r] : r].id : r;
    const double base = (c->qleft - c->y[j]) / (c->qleft - c->anchor) / c->yscale,
                 step = -(c->anchor - c->y[j]) * c->du / c->yscale;
    for (int z = 0; z < nq; ++z)
      c->args[r * nq + z] = fma(step, mid + half * nodes[z], base);
  }
}
/* Implicit balanced range tree over sorted X ranks. Membership is the union
 * of exact compact-X supports, not a Gaussian response-tail approximation.
 * Fully covered nodes emit contiguous donors; empty subtrees never reach the
 * actual kernel-evaluation owner. The original donor identity is retained. */
static void np_cgnn_emit_needed(NPGNNConditionalPrefix *c, int lo, int hi) {
  const int count = c->coverage[hi] - c->coverage[lo];
  if (!count) {
    return;
  }
  if (count == hi - lo) {
    for (int r = lo; r < hi; ++r) {
      int at = c->selectedCount++;
      c->selected[at] = r;
      c->factorIndex[c->order[r].id] = at;
    }
    return;
  }
  const int mid = lo + (hi - lo) / 2;
  np_cgnn_emit_needed(c, lo, mid);
  np_cgnn_emit_needed(c, mid, hi);
}
static void np_cgnn_prepare_tree(NPGNNConditionalPrefix *c) {
  const int n = c->n;
  memset(c->coverage, 0, (size_t)(n + 1) * sizeof(int));
  c->activeCount = 0;
  for (int i = 0; i < n; ++i)
    if (c->active[i]) {
      c->activeIDs[c->activeCount++] = i;
      ++c->coverage[c->first[i]];
      --c->coverage[c->end[i]];
    }
  int overlap = 0, needed = 0;
  for (int r = 0; r < n; ++r) {
    overlap += c->coverage[r];
    c->coverage[r] = needed;
    needed += overlap > 0;
  }
  c->coverage[n] = needed;
  c->selectedCount = 0;
  np_cgnn_emit_needed(c, 0, n);
  if (c->selectedCount != needed)
    c->failed = 1;
}
static double np_cgnn_tree_gap(NPGNNConditionalPrefix *c, int p, int lo, int hi, int z, int nq) {
  const size_t stride = (size_t)(c->selectedCount + 1) * nq;
  const double *v = c->prefix + p * stride, *e = c->correction + p * stride;
  const int a = c->coverage[lo] * nq + z, b = c->coverage[hi] * nq + z;
  return (v[b] - v[a]) + (e[b] - e[a]);
}
static void np_cgnn_tree_rule(NPGNNConditionalPrefix *c, double lo, double hi, int high,
                              double *out) {
  const int nq = high ? 8 : 4, n = c->n, m = c->selectedCount;
#ifdef NP_CF167_TRACE
  np_cgnn_trace_tree_nodes+=(unsigned long)m*nq;
  np_cgnn_trace_pruned_nodes+=(unsigned long)(n-m)*nq;
#endif
  const double *nodes = high ? np_gnn_integral_nodes8 : np_gnn_integral_nodes4,
               *w = high ? np_gnn_integral_weights8 : np_gnn_integral_weights4;
  const double mid = .5 * lo + .5 * hi, half = .5 * hi - .5 * lo;
  if (c->narrow)
    np_cgnn_narrow_args(c, m, nq, 1, 1, mid, half, nodes);
  else
  for (int r = 0; r < m; ++r) {
    const int j = c->order[c->selected[r]].id;
    for (int z = 0; z < nq; ++z)
      c->args[r * nq + z] = (1 + (c->anchor - c->y[j]) * (mid + half * nodes[z])) / c->yscale;
  }
  np_ckernelv(c->ykernel, c->args, m * nq, 0, 0, 1, c->factors, NULL, 0, 0, 1, 1, 0, 0, 0, NULL, NULL);
  for (int a = 0; a < m * nq; ++a)
    c->factors[a] /= c->yscale;

  for (int p = 0; p < c->moments; ++p) {
    double s[8] = {0}, e[8] = {0};
    double *v = c->prefix + (size_t)p * (m + 1) * nq,
           *er = c->correction + (size_t)p * (m + 1) * nq;
    for (int z = 0; z < nq; ++z)
      v[z] = er[z] = 0;
    for (int r = 0; r < m; ++r) {
      const int j = c->order[c->selected[r]].id;
      const double power = c->powers[(size_t)p * n + j];
      for (int z = 0; z < nq; ++z) {
        np_cgnn_compensated(power * c->factors[r * nq + z], s + z, e + z);
        v[(r + 1) * nq + z] = s[z];
        er[(r + 1) * nq + z] = e[z];
      }
    }
  }
  memset(out, 0, (size_t)n * sizeof(double));
  for (int a = 0; a < c->activeCount; ++a) {
    const int i = c->activeIDs[a];
    for (int z = 0; z < nq; ++z) {
      double fit = c->coef[c->moments * i] * np_cgnn_tree_gap(c, 0, c->first[i], c->end[i], z, nq);
      for (int p = 1; p < c->moments; ++p)
        fit = fma(c->coef[c->moments * i + p],
                  np_cgnn_tree_gap(c, p, c->first[i], c->end[i], z, nq), fit);
      fit -= c->self[i] * c->factors[c->factorIndex[i] * nq + z];
      out[i] += half * w[z] * fit * fit;
    }
  }
}
static double np_cgnn_prefixgap(NPGNNConditionalPrefix *c, int p, int lo, int hi, int z, int nq) {
  const size_t stride = (size_t)(c->n + 1) * nq;
  const double *v = c->prefix + p * stride, *e = c->correction + p * stride;
  return (v[hi * nq + z] - v[lo * nq + z]) + (e[hi * nq + z] - e[lo * nq + z]);
}
static void np_cgnn_rule(NPGNNConditionalPrefix *c, double lo, double hi, int high, double *out) {
  if (c->stable || c->local.degree) {
    np_cgnn_local_rule(c, lo, hi, high, 0, out);
    return;
  }
  /* Exact root certificate, not a density/cost heuristic: if the support union
   * contains every donor, the existing dense arithmetic is already minimal.
   * The tree was traversed and remains enabled. No statistical method changes.
   * Otherwise every tree-pruned donor is omitted by the compacted owner. */
  if (c->tree && c->selectedCount != c->n) {
    np_cgnn_tree_rule(c, lo, hi, high, out);
    return;
  }
  const int nq = high ? 8 : 4, n = c->n;
  const double *nodes = high ? np_gnn_integral_nodes8 : np_gnn_integral_nodes4;
  const double *w = high ? np_gnn_integral_weights8 : np_gnn_integral_weights4;
  const double mid = .5 * lo + .5 * hi, half = .5 * hi - .5 * lo;
  if (c->narrow)
    np_cgnn_narrow_args(c, n, nq, 0, 0, mid, half, nodes);
  else
  for (int j = 0; j < n; ++j)
    for (int z = 0; z < nq; ++z)
      c->args[j * nq + z] = (1 + (c->anchor - c->y[j]) * (mid + half * nodes[z])) / c->yscale;
  np_ckernelv(c->ykernel, c->args, n * nq, 0, 0, 1, c->factors, NULL, 0, 0, 1, 1, 0, 0, 0, NULL, NULL);
  for (int at = 0; at < n * nq; ++at)
    c->factors[at] /= c->yscale;
  for (int p = 0; p < c->moments; ++p) {
    double s[8] = {0}, e[8] = {0};
    double *v = c->prefix + (size_t)p * (n + 1) * nq,
           *er = c->correction + (size_t)p * (n + 1) * nq;
    for (int z = 0; z < nq; ++z)
      v[z] = er[z] = 0;
    for (int r = 0; r < n; ++r) {
      const int j = c->order[r].id;
      const double power = c->powers[(size_t)p * n + j];
      for (int z = 0; z < nq; ++z) {
        np_cgnn_compensated(power * c->factors[j * nq + z], s + z, e + z);
        v[(r + 1) * nq + z] = s[z];
        er[(r + 1) * nq + z] = e[z];
      }
    }
  }

  for (int i = 0; i < n; ++i) {
    out[i] = 0;
    if (!c->active[i])
      continue;
    for (int z = 0; z < nq; ++z) {
      double fit = c->coef[c->moments * i] * np_cgnn_prefixgap(c, 0, c->first[i], c->end[i], z, nq);
      for (int p = 1; p < c->moments; ++p)
        fit = fma(c->coef[c->moments * i + p],
                  np_cgnn_prefixgap(c, p, c->first[i], c->end[i], z, nq), fit);
      fit -= c->self[i] * c->factors[i * nq + z];
      out[i] += half * w[z] * fit * fit;
    }
  }
}
static void np_cgnn_gaussian_derivative_constants(int kernel, double logG[17]) {
 double lower[23]={0},upper[23]={0};
 if(kernel==0)lower[0]=1;
 if(kernel==1){lower[0]=1.5;lower[2]=-.5;}
 if(kernel==2){lower[0]=1.875;lower[2]=-1.25;lower[4]=.125;}
 if(kernel==3){lower[0]=2.1875;lower[2]=-2.1875;
  lower[4]=.4375;lower[6]=-.02083333333;}
 for(int p=0;p<23;++p)upper[p]=lower[p];
 for(int r=0;r<=16;++r){
  double sum=0;
  for(int p=0;p<=6+r;++p){
   const double a=fmax(fabs(lower[p]),fabs(upper[p]));
   const double monomial=p?exp(.5*p*(log(2.0*p)-1)):1;
   sum+=a*monomial;
  }
  /* Coefficient operations below are outward enclosed. This small inflation
   * also covers ordinary rounding in this finite positive sum and libm;
   * it is not a change to the quadrature accuracy target. */
  logG[r]=log(sum)-.5*log(2*M_PI)+log1p(256*DBL_EPSILON);
  if(r==16)break;
  double nextLower[23]={0},nextUpper[23]={0};
  for(int p=0;p<=7+r;++p){
   const double a=p+1<23?lower[p+1]:0,b=p+1<23?upper[p+1]:0;
   const double c=p?lower[p-1]:0,d=p?upper[p-1]:0;
   nextLower[p]=nextafter(nextafter((p+1)*a,-INFINITY)-d,-INFINITY);
   nextUpper[p]=nextafter(nextafter((p+1)*b,INFINITY)-c,INFINITY);
  }
  for(int p=0;p<23;++p){lower[p]=nextLower[p];upper[p]=nextUpper[p];}
 }
}
/* The reciprocal response is K((1+(anchor-y)*u)/yscale)/yscale.
 * Signed canonical influence rows multiply the squared-fit bound by maxL1^2.
 * Retain exact whole-support integration and the existing accuracy targets. */
static double np_cgnn_gaussian_interval_bound(NPGNNConditionalPrefix *c,
                                             double lo, double hi) {
 const double *y=c->y,*logG=c->logG,anchor=c->anchor;
 const int n=c->n;
 const double width=nextafter(hi-lo,INFINITY),maxU=fmax(fabs(lo),fabs(hi));
 const double logWidth=log(width);
 double maxima[17];for(int r=0;r<=16;++r)maxima[r]=-INFINITY;
 for(int j=0;j<n;++j){
  const double slope=anchor-y[j];
  const double magnitude=fabs(anchor)+fabs(y[j]);
  double slopeUpper=nextafter((fabs(slope)+DBL_EPSILON*magnitude)/c->yscale,INFINITY);
  double a=fma(slope,lo,1)/c->yscale,b=fma(slope,hi,1)/c->yscale;
  double roundoff=8*DBL_EPSILON*(1+magnitude*maxU)/c->yscale;
  if(c->narrow) {
    const double base=(c->qleft-y[j])/(c->qleft-anchor)/c->yscale,
                 step=-slope*c->du/c->yscale;
    slopeUpper=nextafter((fabs(slope)+DBL_EPSILON*magnitude)*c->du/c->yscale,INFINITY);
    a=fma(step,lo,base);b=fma(step,hi,base);
    roundoff=16*DBL_EPSILON*(fabs(base)+slopeUpper*maxU);
  }
  const double nearest=(a<=0&&b>=0)||(a>=0&&b<=0)?0:
   fmax(0,fmin(fabs(a),fabs(b))-roundoff);
  const double decay=-.25*nearest*nearest;
  maxima[0]=fmax(maxima[0],decay);
  if(slopeUpper>0){
   const double step=log(slopeUpper)+logWidth;
   double term=decay;
   for(int r=1;r<=16;++r){term+=step;maxima[r]=fmax(maxima[r],term);}
  }
 }
 double terms[17],largest=-INFINITY,choose=1;
 for(int r=0;r<=16;++r){
  terms[r]=log(choose)+logG[r]+logG[16-r]+maxima[r]+maxima[16-r];
  largest=fmax(largest,terms[r]);choose*=((double)(16-r))/(r+1);
 }
 if(largest==-INFINITY)return 0;
 if(!R_FINITE(largest))return INFINITY;
 double sum=0;for(int r=0;r<=16;++r)sum+=exp(terms[r]-largest);
 const double logBound=2*(log(c->maxL1)-log(c->yscale))+logWidth+4*lgamma(9.0)-log(17.0)-3*lgamma(17.0)
  +largest+log(sum)+log1p(1024*DBL_EPSILON);
 return logBound<log(DBL_MIN)?DBL_MIN:nextafter(exp(logBound),INFINITY);
}

static double np_cgnn_missed_peak_bound(NPGNNConditionalPrefix *c, double lo, double hi) {
  if (c->ykernel != 0 || c->narrow)
    return np_cgnn_gaussian_interval_bound(c, lo, hi);
  const double width = hi - lo,
               variation = fmax(fabs(c->anchor - c->minY), fabs(c->anchor - c->maxY)) / c->yscale;
  if (variation * width <= 4)
    return 0;
  double logMax = R_NegInf;
  for (int j = 0; j < c->n; ++j) {
    const double slope = (c->anchor - c->y[j]) / c->yscale;
    if (fabs(slope) * width <= 4)
      continue;
    const double a = 1 / c->yscale + slope * lo, b = 1 / c->yscale + slope * hi;
    const double nearest = (a <= 0 && b >= 0) || (b <= 0 && a >= 0) ? 0 : fmin(fabs(a), fabs(b));
    logMax = fmax(logMax, (log(allck[0](0)) - log(c->yscale)) - .5 * nearest * nearest);
  }
  if (logMax == R_NegInf)
    return 0;
  /* Signed LP rows need their actual L1 bound, not the LC value one. */
  const double logFull = (log(allck[0](0)) - log(c->yscale));
  const double bound = 2 * log(c->maxL1) + log(width) + log(4.0) + logFull + logMax +
                       log1p(.5 * exp(logMax - logFull));
  return bound < log(DBL_MIN) ? DBL_MIN : exp(bound);
}
static void np_cgnn_integrate(NPGNNConditionalPrefix *c, double lo, double hi, int depth,
                              double *out) {
  if (c->failed)
    return;
  ++c->visits;
  if(++c->interval_visits > c->work_limit){
    int row=0;while(row<c->n && !c->active[row])row++;
    np_conditional_failure_record(NP_CONDITIONAL_WORK_EXHAUSTED,row+1);
    c->failed=1;return;
  }
  const int nf = c->folds;
  double *low = c->recursive + (size_t)depth * 3 * nf, *high = low + nf, *other = high + nf;
  /* A certified local remainder may spend the existing relative budget
   * only against a lower bound on this piece's true mean-fold integral. */
  const int certificate_domain=!c->narrow && c->ykernel == 0 && c->yscale == 1.0 &&
    np_lp_engine_extern == NP_LP_ENGINE_SCALAR && KERNEL_reg_extern == 4;
  double bound=R_PosInf;
  if(certificate_domain) {
    const double log_bound=np_cgnn_gauss4_log_bound(lo,hi,c->anchor-c->maxY,
      c->anchor-c->minY,(double)c->activeCount/nf*c->maxL1*c->maxL1);
    bound=exp(log_bound);
    /* Preserve the established absolute-only shortcut without extra work. */
    if(log_bound <= log(ldexp(1e-12/(8*c->n),-depth))) {
#ifdef NP_CF167_TRACE
      ++np_cgnn_trace_certificates;
#endif
      np_cgnn_rule(c,lo,hi,0,out);
      for(int i=0;i<nf;++i)if(!R_FINITE(out[i])){c->failed=1;return;}
      c->acceptedError+=bound;
      return;
    }
  }
  np_cgnn_rule(c, lo, hi, 0, low);
  if(certificate_domain) {
    long double sum=0;
    for(int i=0;i<nf;++i) {
      if(!R_FINITE(low[i])){c->failed=1;return;}
      sum+=(long double)low[i]/nf;
    }
    const double accumulation=(2.0*nf+8.0)*DBL_EPSILON;
    const double lower=fmax(0.0,nextafter((double)sum/(1.0+accumulation),0.0)-bound);
    if(bound <= ldexp(1e-12/(8*c->n),-depth)+1e-10*lower) {
#ifdef NP_CF167_TRACE
      ++np_cgnn_trace_certificates;
#endif
      memcpy(out,low,(size_t)nf*sizeof(double));
      c->acceptedError+=bound;
      return;
    }
  }
  np_cgnn_rule(c, lo, hi, 1, high);

  const double budget = ldexp(1e-12 / (8 * c->n), -depth) / (c->narrow ? c->du : 1.0),
               unseen = np_cgnn_missed_peak_bound(c, lo, hi);
  int pass = R_FINITE(unseen) && unseen <= budget / 4;
  double err = 0, value = 0;
  for (int i = 0; i < nf; ++i) {
    if (!R_FINITE(high[i]) || !R_FINITE(low[i])) {
      c->failed = 1;
      return;
    }
    if (!c->active[i])
      continue;
    const double e = fabs(low[i] - high[i]) + unseen;
    err += e / nf;
    value += fabs(high[i]) / nf;
  }
  /* I1 is the mean over deleted folds. Sum ABSOLUTE error estimates before
   * comparing at that same scale, so errors cannot cancel. Retain the same
   * global absolute/relative targets and the independent missed-peak gate. */
  if (!R_FINITE(err) || !R_FINITE(value) || !R_FINITE(unseen)) {
    c->failed = 1;
    return;
  }
  if (err > budget + 1e-10 * value)
    pass = 0;
  if (pass) {
    for (int i = 0; i < nf; ++i)
      out[i] = high[i];
    c->acceptedError += err;
    return;
  }
  const double mid = .5 * lo + .5 * hi;
  if (depth >= 20 || !(lo < mid && mid < hi)) {
    int row=0;while(row<c->n && !c->active[row])row++;
    np_conditional_failure_record(NP_CONDITIONAL_WORK_EXHAUSTED,row+1);
    c->failed=1;return;
  }

  np_cgnn_integrate(c, lo, mid, depth + 1, out);
  if (c->failed)
    return;
  np_cgnn_integrate(c, mid, hi, depth + 1, other);
  if (c->failed)
    return;
  for (int i = 0; i < nf; ++i)
    out[i] += other[i];
}

/* Compact responses are polynomial on reciprocal-radius intervals split at
 * every support crossing. A degree-d response needs d+1 Gauss nodes for its
 * square. This kernel-family representation has no refinement/failure switch;
 * the native kernel remains the source of all floating-point factors. */
static const double np_cgnn_nodes1[] = {0}, np_cgnn_weights1[] = {2};
static const double np_cgnn_nodes3[] = {-.7745966692414833770,0,.7745966692414833770};
static const double np_cgnn_weights3[] = {.5555555555555555556,.8888888888888888889,.5555555555555555556};
static const double np_cgnn_nodes5[] = {-.9061798459386639928,-.5384693101056830910,0,.5384693101056830910,.9061798459386639928};
static const double np_cgnn_weights5[] = {.2369268850561890875,.4786286704993664680,.5688888888888888889,.4786286704993664680,.2369268850561890875};
static const double np_cgnn_nodes7[] = {-.9491079123427585245,-.7415311855993944399,-.4058451513773971669,0,.4058451513773971669,.7415311855993944399,.9491079123427585245};
static const double np_cgnn_weights7[] = {.1294849661688696933,.2797053914892766679,.3818300505051189450,.4179591836734693878,.3818300505051189450,.2797053914892766679,.1294849661688696933};
static const double np_cgnn_nodes9[] = {-.9681602395076260898,-.8360311073266357943,-.6133714327005903973,-.3242534234038089290,0,.3242534234038089290,.6133714327005903973,.8360311073266357943,.9681602395076260898};
static const double np_cgnn_weights9[] = {.08127438836157441197,.1806481606948574041,.2606106964029354623,.3123470770400028401,.3302393550012597632,.3123470770400028401,.2606106964029354623,.1806481606948574041,.08127438836157441197};

static void np_cgnn_local_rule(NPGNNConditionalPrefix *c, double lo, double hi,
                                int high, int compact, double *out) {
  const int nq = compact ? (c->ykernel == 8 ? 1 : 2*(c->ykernel-3)+1) : (high ? 8 : 4);
  const int pruned = c->tree && c->selectedCount != c->n;
  const int m = pruned ? c->selectedCount : c->n;
  const double *nodes = compact ? (nq == 1 ? np_cgnn_nodes1 : nq == 3 ? np_cgnn_nodes3 :
    nq == 5 ? np_cgnn_nodes5 : nq == 7 ? np_cgnn_nodes7 : np_cgnn_nodes9) :
    (high ? np_gnn_integral_nodes8 : np_gnn_integral_nodes4);
  const double *weights = compact ? (nq == 1 ? np_cgnn_weights1 : nq == 3 ? np_cgnn_weights3 :
    nq == 5 ? np_cgnn_weights5 : nq == 7 ? np_cgnn_weights7 : np_cgnn_weights9) :
    (high ? np_gnn_integral_weights8 : np_gnn_integral_weights4);
  const double mid = .5*lo+.5*hi, half = .5*hi-.5*lo;
  if (c->narrow)
    np_cgnn_narrow_args(c, m, nq, 1, pruned, mid, half, nodes);
  else
    for (int r = 0; r < m; ++r) {
      const int j = c->order[pruned ? c->selected[r] : r].id;
      for (int z = 0; z < nq; ++z)
        c->args[r*nq+z] = (1+(c->anchor-c->y[j])*(mid+half*nodes[z]))/c->yscale;
    }
  np_ckernelv(c->ykernel,c->args,m*nq,0,0,1,c->factors,NULL,0,0,1,1,0,0,0,NULL,NULL);
  int retained = 0;
  for (int r = 0; r < m; ++r) {
    int nonzero = 0;
    for (int z = 0; z < nq; ++z) {
      c->factors[r*nq+z] /= c->yscale;
      nonzero |= c->factors[r*nq+z] != 0.0;
    }
    if (nonzero) {
      c->local.ranks[retained] = pruned ? c->selected[r] : r;
      c->local.slots[retained] = r;
      if(c->stable && c->full_rows){
        c->stable_moments.ranks[retained]=c->local.ranks[retained];
        c->stable_moments.slots[retained]=r;
      }
      retained++;
    }
  }
#ifdef NP_CF167_TRACE
  np_cgnn_trace_tree_nodes += (unsigned long)m*nq;
  np_cgnn_trace_pruned_nodes += (unsigned long)(c->n-m)*nq;
#endif
  if(c->stable && c->full_rows)
    np_cgnn_stable_evaluate(&c->stable_moments,c->ones,c->factors,nq,retained,c->active,c->radius);
  if(!c->stable || c->deficient_rows)
    np_cgnn_local_evaluate(&c->local,c->moments,c->powers,c->coef,c->factors,
                          nq,retained,c->active);
  memset(out,0,(size_t)c->n*sizeof(double));
  for (int i = 0; i < c->n; ++i) if (c->active[i])
    for (int z = 0; z < nq; ++z) {
      const double value = c->stable && c->radius[i]>0 ?
        c->stable_moments.fit[(size_t)i*nq+z] : c->local.fit[(size_t)i*nq+z];
      out[i] += half*weights[z]*value*value;
    }
}

static void np_cgnn_compact_rule(NPGNNConditionalPrefix *c, double lo, double hi, double *out) {
  if (c->stable || c->local.degree) {
    np_cgnn_local_rule(c, lo, hi, 1, 1, out);
    return;
  }
  const int nq = c->ykernel == 8 ? 1 : 2 * (c->ykernel - 3) + 1, n = c->n;
  const int pruned = c->tree && c->selectedCount != n;
  const int m = pruned ? c->selectedCount : n;
  const double *nodes = nq == 1 ? np_cgnn_nodes1 : nq == 3 ? np_cgnn_nodes3 :
      nq == 5 ? np_cgnn_nodes5 : nq == 7 ? np_cgnn_nodes7 : np_cgnn_nodes9;
  const double *weights = nq == 1 ? np_cgnn_weights1 : nq == 3 ? np_cgnn_weights3 :
      nq == 5 ? np_cgnn_weights5 : nq == 7 ? np_cgnn_weights7 : np_cgnn_weights9;
  const double mid = .5 * lo + .5 * hi, half = .5 * hi - .5 * lo;
  if (c->narrow)
    np_cgnn_narrow_args(c, m, nq, 1, pruned, mid, half, nodes);
  else
  for (int r = 0; r < m; ++r) {
    const int j = c->order[pruned ? c->selected[r] : r].id;
    for (int z = 0; z < nq; ++z)
      c->args[r * nq + z] = (1 + (c->anchor - c->y[j]) * (mid + half * nodes[z])) / c->yscale;
  }
  np_ckernelv(c->ykernel, c->args, m * nq, 0, 0, 1, c->factors, NULL, 0, 0, 1, 1, 0, 0, 0, NULL, NULL);
  int retained = 0;
  c->responseRank[0] = 0;
  for (int r = 0; r < m; ++r) {
    int nonzero = 0;
    for (int z = 0; z < nq; ++z) {
      c->factors[r * nq + z] /= c->yscale;
      nonzero |= c->factors[r * nq + z] != 0.0;
    }
    c->responseActive[r] = nonzero;
    if (nonzero)
      c->responseSlots[retained++] = r;
    c->responseRank[r + 1] = retained;
  }
  /* responseRank maps selected-X slots to response-prefix ranks, never donor
   * identities. factorIndex/orderRank still address uncompressed factors. */
  const size_t stride = (size_t)(retained + 1) * nq;
  for (int p = 0; p < c->moments; ++p) {
    double *v = c->prefix + p * stride, *e = c->correction + p * stride;
    for (int z = 0; z < nq; ++z)
      v[z] = e[z] = 0.0;
    double sum[9] = {0}, error[9] = {0};
    for (int r = 0; r < retained; ++r) {
      const int slot = c->responseSlots[r];
      const int j = c->order[pruned ? c->selected[slot] : slot].id;
      const double power = c->powers[(size_t)p * n + j];
      for (int z = 0; z < nq; ++z) {
        np_cgnn_compensated(power * c->factors[slot * nq + z], sum + z, error + z);
        v[(r + 1) * nq + z] = sum[z];
        e[(r + 1) * nq + z] = error[z];
      }
    }
  }
  memset(out, 0, (size_t)n * sizeof(double));
  for (int a = 0; a < (c->tree ? c->activeCount : n); ++a) {
    const int i = c->tree ? c->activeIDs[a] : a;
    if (!c->active[i])
      continue;
    const int first = c->responseRank[pruned ? c->coverage[c->first[i]] : c->first[i]],
              end = c->responseRank[pruned ? c->coverage[c->end[i]] : c->end[i]],
              self = pruned ? c->factorIndex[i] : c->orderRank[i];
    if (first == end && !c->responseActive[self])
      continue; /* All contributing computed factors are exactly zero. */
    for (int z = 0; z < nq; ++z) {
      double value = 0.0;
      for (int p = 0; p < c->moments; ++p) {
        const double *v = c->prefix + p * stride, *e = c->correction + p * stride;
        const double gap = (v[end * nq + z] - v[first * nq + z]) +
                           (e[end * nq + z] - e[first * nq + z]);
        value = fma(c->coef[(size_t)c->moments * i + p], gap, value);
      }
      value -= c->self[i] * c->factors[self * nq + z];
      out[i] += half * weights[z] * value * value;
    }
  }
}

static int np_cgnn_compare_cut(const void *a, const void *b) {
  const double x = *(const double *)a, y = *(const double *)b;
  return (x > y) - (x < y);
}
static void np_cgnn_compact_integrate(NPGNNConditionalPrefix *c, double lo, double hi, double *out) {
  const double support = c->ykernel == 8 ? 1.0 : sqrt(5.0);
  const int pruned = c->tree && c->selectedCount != c->n;
  const int m = pruned ? c->selectedCount : c->n;
  int count = 0;
  c->cuts[count++] = lo;
  c->cuts[count++] = hi;
  for (int r = 0; r < m; ++r) {
    const int j = c->order[pruned ? c->selected[r] : r].id;
    const double slope = c->anchor - c->y[j];
    if (slope == 0.0)
      continue;
    for (int side = -1; side <= 1; side += 2) {
      double u = (side * support * c->yscale - 1.0) / slope;
      if (c->narrow) {
        const double base = (c->qleft - c->y[j]) / (c->qleft - c->anchor),
                     step = -slope * c->du;
        u = (side * support * c->yscale - base) / step;
      }
      if (lo < u && u < hi)
        c->cuts[count++] = u;
    }
  }
  qsort(c->cuts, count, sizeof(double), np_cgnn_compare_cut);
  double *one = c->recursive, *error = one + c->n;
  memset(out, 0, (size_t)c->n * sizeof(double));
  memset(error, 0, (size_t)c->n * sizeof(double));
  for (int at = 1; at < count; ++at) {
    if (!(c->cuts[at - 1] < c->cuts[at]))
      continue;
    np_progress_bandwidth_loop_step();
    np_cgnn_compact_rule(c, c->cuts[at - 1], c->cuts[at], one);
    for (int i = 0; i < c->n; ++i)
      np_cgnn_compensated(one[i], out + i, error + i);
  }
  for (int i = 0; i < c->n; ++i) {
    out[i] += error[i];
    if (!R_FINITE(out[i]))
      c->failed = 1;
  }
}

/* For unit-scale uniform Y, the anchor donor is identically zero at the
 * kernel's open support boundary. Every other donor has bounded reciprocal
 * support between 0 and -2/(anchor-y). Their union is an exact finite domain,
 * even when a response endpoint has zero radius. This is not a tail cutoff. */
static int np_cgnn_uniform_finite_panel(NPGNNConditionalPrefix *c, double *lo, double *hi) {
  if (c->ykernel != 8 || c->yscale != 1.0)
    return -1;
  double left = 0.0, right = 0.0;
  const int pruned = c->tree && c->selectedCount != c->n;
  const int m = pruned ? c->selectedCount : c->n;
  for (int r = 0; r < m; ++r) {
    const int j = c->order[pruned ? c->selected[r] : r].id;
    const double slope = c->anchor - c->y[j];
    if (slope == 0.0)
      continue;
    const double edge = -2.0 / slope;
    if (!R_FINITE(edge))
      return -1;
    left = fmin(left, edge);
    right = fmax(right, edge);
  }
  *lo = fmax(*lo, left);
  *hi = fmin(*hi, right);
  return *lo < *hi;
}

/* A mathematically exact representation boundary, not a failure fallback. */
static int np_cgnn_prefix_admitted(void) {
  return BANDWIDTH_den_extern == BW_GEN_NN && num_reg_continuous_extern == 1 &&
         num_var_continuous_extern == 1 && num_reg_unordered_extern == 0 &&
         num_reg_ordered_extern == 0 && num_var_unordered_extern == 0 &&
         num_var_ordered_extern == 0 &&
         (KERNEL_reg_extern >= 4 && KERNEL_reg_extern <= 8) && KERNEL_den_extern >= 0 && KERNEL_den_extern <= 8 &&
         int_cxker_bound_extern == 0 && int_cyker_bound_extern == 0;
}

/* Only completed objectives contribute notices. Never signal in a callback
 * or rank-owned loop: warn=2 must not strand ranks at later collectives. */
static double np_cgnn_notice_caps, np_cgnn_notice_error;
void np_conditional_gnn_notice_reset(void) {
  np_cgnn_notice_caps = np_cgnn_notice_error = 0.0;
#ifdef NP_CF167_TRACE
  np_cgnn_trace_tree_nodes=np_cgnn_trace_pruned_nodes=0;
  np_cgnn_trace_projected_nodes=np_cgnn_trace_certificates=0;
#endif
}
void np_conditional_gnn_notice_emit(void) {
  const double caps = np_cgnn_notice_caps, estimate = np_cgnn_notice_error;
#ifdef NP_CF167_TRACE
  REprintf("CF167 owner rank=%d tree_nodes=%lu pruned=%lu projected_nodes=%lu certificates=%lu\n",
    np_cgnn_trace_rank(),np_cgnn_trace_tree_nodes,np_cgnn_trace_pruned_nodes,
    np_cgnn_trace_projected_nodes,np_cgnn_trace_certificates);
#endif
  np_conditional_gnn_notice_reset();
#ifdef MPI2
  if (my_rank != 0 && !np_mpi_local_regression_active())
    return;
#endif
  if (caps > 0.0)
    warning("conditional GNN cv.ls completed the full integral, but %.0f refinement pieces reached "
            "the quadrature limit; maximum estimated absolute I1 error %.3g (targets: absolute "
            "1e-12, relative 1e-10). Inspect bandwidth sensitivity before interpreting small "
            "objective differences",
            caps, estimate);
}

typedef struct {
  NPGNNConditionalPrefix integral;
  NPConditionalQRDeleted qr;
  NPConditionalXRowCtx xctx;
  NPConditionalYRowCtx yctx;
  NPGNNIntegralGeometry geometry;
  double *vsf, *y, *row, *yrow, *folds, *pieces, *out;
  double score;
  int status;
} NPGNNConditionalCall;

static void np_cgnn_cleanup(void *raw, Rboolean jump) {
  NPGNNConditionalCall *a = (NPGNNConditionalCall *)raw;
  NPGNNConditionalPrefix *c = &a->integral;
  (void)jump;
  np_cqr_deleted_clear(&a->qr);
  np_cgnn_stable_clear(&c->stable_moments);
  free(c->ones);free(c->coefficient_work);
  np_conditional_xrow_ctx_clear(&a->xctx);
  np_conditional_yrow_ctx_clear(&a->yctx);
  np_gnn_integral_geometry_clear(&a->geometry);
  free(a->y);
  free(a->row);
  free(a->yrow);
  free(a->folds);
  free(a->pieces);
  free(a->out);
  free(c->active);
  free(c->first);
  free(c->end);
  free(c->orderRank);
  free(c->x);
  free(c->order);
  free(c->coef);
  free(c->self);
  free(c->powers);
  free(c->args);
  free(c->factors);
  free(c->prefix);
  free(c->correction);
  free(c->recursive);
  free(c->coverage);
  free(c->selected);
  free(c->factorIndex);
  free(c->activeIDs);
  free(c->responseRank);
  free(c->responseActive);
  free(c->responseSlots);
  free(c->cuts);
  np_cgnn_local_clear(&c->local);
  np_glp_cv_clear_extern();
}

static void *np_cgnn_calloc(size_t n, size_t width) {
  size_t bytes;
  return np_size_mul_checked(n, width, &bytes) ? calloc(1, bytes) : NULL;
}

static SEXP np_cgnn_body(void *raw) {
  NPGNNConditionalCall *a = (NPGNNConditionalCall *)raw;
  NPGNNConditionalPrefix *c = &a->integral;
  NPConditionalXRowCtx *x = &a->xctx;
  NPGNNIntegralGeometry *g = &a->geometry;
  const int n = num_obs_train_extern;
  const int compact = KERNEL_den_extern >= 4;
  const int maxnq = compact ? 9 : 8;
  const int uniform = KERNEL_reg_extern == 8;
  const int foldPlanes = uniform ? 4 : 2;
  int fail = 0, first = 0, count = n;
#ifdef MPI2
  const int parallel = np_objective_outer_rows_enabled(1);
#endif
  a->status = 1;
  if (n < 3 || n > INT_MAX / maxnq)
    return R_NilValue;
  fail = np_conditional_xrow_ctx_prepare_ctx(a->vsf, &np_conditional_deleted_identity_geometry, x);
  if (!fail)
    fail = np_conditional_yrow_ctx_prepare_ctx(a->vsf, OP_NORMAL,
                                               &np_conditional_deleted_identity_geometry, &a->yctx);
  const int scalar = np_lp_engine_extern == NP_LP_ENGINE_SCALAR;
  const int terms = scalar ? 1 : np_glp_cv_cache.nterms;
  if (terms < 1 || (size_t)terms > INT_MAX / (size_t)n)
    fail = 1;
  c->n = c->folds = n;
  c->moments = fail ? 0 : terms;
  c->stable=!scalar && terms>1 && (KERNEL_reg_extern==4 || uniform);
  c->coefficient_count=c->moments*n;
  if(c->stable){
    if((size_t)n>INT_MAX/((size_t)3*terms+1))fail=1;
    else c->coefficient_count=n*(3*terms+1);
  }
  /* Allocate 256 recursive visits to each geometry interval up front.
   * Both deletion branches and tails spend that same interval allocation.
   * Total work is bounded by 256*geometry.count, independently of MPI ranks. */
  c->work_limit=256;
  c->ykernel = KERNEL_den_extern;
  if (c->ykernel > 0 && c->ykernel < 4)
    np_cgnn_gaussian_derivative_constants(c->ykernel, c->logG);
  c->maxL1 = 1.0;
  c->tree = int_TREE_X == NP_TREE_TRUE;
  size_t prefix_count = 0, recursive_count = 0;
  /* Check products before forming them, including on 32-bit size_t hosts. */
  if (!fail && (!np_size_mul_checked(c->moments, (size_t)n + 1, &prefix_count) ||
                !np_size_mul_checked(prefix_count, maxnq, &prefix_count) ||
                !np_size_mul_checked(n, compact ? 2 : 63, &recursive_count)))
    fail = 1;
#define NP_CGNN_ALLOC(ptr, length)                                                                 \
  do {                                                                                             \
    if (!fail) {                                                                                   \
      (ptr) = np_cgnn_calloc((length), sizeof(*(ptr)));                                            \
      if (!(ptr))                                                                                  \
        fail = 1;                                                                                  \
    }                                                                                              \
  } while (0)
  NP_CGNN_ALLOC(a->y, n);
  NP_CGNN_ALLOC(a->row, n);
  NP_CGNN_ALLOC(a->yrow, n);
  NP_CGNN_ALLOC(a->out, n);
  NP_CGNN_ALLOC(a->folds, (size_t)foldPlanes * n);
  NP_CGNN_ALLOC(c->active, n);
  NP_CGNN_ALLOC(c->first, n);
  NP_CGNN_ALLOC(c->end, n);
  NP_CGNN_ALLOC(c->x, n);
  NP_CGNN_ALLOC(c->order, n);
  if (uniform || compact) {
    NP_CGNN_ALLOC(c->orderRank, n);
  }
  NP_CGNN_ALLOC(c->coef, c->coefficient_count);
  if(c->stable && !fail){
    c->stable_coef=(NPConditionalQRAcc *)(c->coef+(size_t)terms*n);
    c->radius=c->coef+(size_t)3*terms*n;
    NP_CGNN_ALLOC(c->coefficient_work,terms);NP_CGNN_ALLOC(c->ones,n);
    if(!fail)for(int j=0;j<n;j++)c->ones[j]=1.;
  }
  NP_CGNN_ALLOC(c->self, n);
  NP_CGNN_ALLOC(c->powers, (size_t)c->moments * n);
  NP_CGNN_ALLOC(c->args, (size_t)maxnq * n);
  NP_CGNN_ALLOC(c->factors, (size_t)maxnq * n);
  NP_CGNN_ALLOC(c->prefix, prefix_count);
  NP_CGNN_ALLOC(c->correction, prefix_count);
  NP_CGNN_ALLOC(c->recursive, recursive_count);
  if (compact) {
    NP_CGNN_ALLOC(c->responseRank, n + 1);
    NP_CGNN_ALLOC(c->responseActive, n);
    NP_CGNN_ALLOC(c->responseSlots, n);
    NP_CGNN_ALLOC(c->cuts, (size_t)2 * n + 2);
  }
  if (c->tree) {
    NP_CGNN_ALLOC(c->coverage, n + 1);
    NP_CGNN_ALLOC(c->selected, n);
    NP_CGNN_ALLOC(c->factorIndex, n);
    NP_CGNN_ALLOC(c->activeIDs, n);
  }
  if (!fail) {
    for (int j = 0; j < n; ++j)
      a->y[int_TREE_Y == NP_TREE_TRUE ? ipt_extern_Y[j] : j] =
          matrix_Y_continuous_train_extern[0][j];
    fail = np_gnn_integral_geometry_prepare(g, a->y, n, a->yctx.vsfy[0], 1);
  }
  if (!fail && g->count > INT_MAX / 3)
    fail = 1;
  NP_CGNN_ALLOC(a->pieces, (size_t)3 * g->count);
#undef NP_CGNN_ALLOC
#ifdef MPI2
  if (np_conditional_outer_preflight_failed(parallel, fail))
    return R_NilValue;
  np_objective_outer_owned_rows(0, n, parallel, &first, &count);
#else
  if (fail)
    return R_NilValue;
#endif
  c->y = a->y;
  c->yscale = g->scale;
  c->minY = g->sorted[0];
  c->maxY = g->sorted[n - 1];
  double xmin = matrix_X_continuous_train_extern[0][0], xmax = xmin;
  for (int j = 1; j < n; ++j) {
    xmin = fmin(xmin, matrix_X_continuous_train_extern[0][j]);
    xmax = fmax(xmax, matrix_X_continuous_train_extern[0][j]);
  }
  const double center = .5 * xmin + .5 * xmax, scale = .5 * xmax - .5 * xmin;
  if (!(scale > 0.0) || !R_FINITE(scale))
    fail = 1;
  for (int j = 0; j < n && !fail; ++j) {
    const int id = int_TREE_X == NP_TREE_TRUE ? ipt_extern_X[j] : j;
    const double z = (matrix_X_continuous_train_extern[0][j] - center) / scale;
    c->x[id] = uniform ? z : matrix_X_continuous_train_extern[0][j];
    c->order[j] = (NPGNNConditionalOrder){
        matrix_X_continuous_train_extern[0][j], id};
    for (int t = 0; t < terms; ++t) {
      const double b = scalar ? 1.0 : x->basis[t][j];
      c->powers[(size_t)t * n + id] = b;
    }
  }
  qsort(c->order, n, sizeof(*c->order), np_cgnn_compare);
  if (uniform || compact)
    for (int r = 0; r < n; ++r)
      c->orderRank[c->order[r].id] = r;
  const double k0 = allck[KERNEL_reg_extern](0.0);
  for (int i = first; i < first + count && !fail; ++i) {
    np_progress_bandwidth_loop_step();
    if ((c->stable ? np_conditional_deleted_from_ctx_core(x,&a->qr,i,1,0,a->row,NULL) :
        np_conditional_xrow_from_ctx(x, i, a->row)) ||
        (!scalar && !np_cgnn_deleted_support_sufficient(
          c->order, x->kw, int_TREE_X == NP_TREE_TRUE ? ipt_lookup_extern_X : NULL,
          n, i, terms)) ||
        np_conditional_yrow_from_ctx(&a->yctx, i, a->yrow)) {
      fail = 1;
      break;
    }
    const int pos = int_TREE_X == NP_TREE_TRUE ? ipt_lookup_extern_X[i] : i;
    const double h = x->matrix_bandwidth_eval_one[0][0];
    if (uniform) {
      /* Reuse the canonical strict-distance weight mask. Transformed endpoint
       * comparisons can change membership at the uniform kernel's jump.
       * Include the query here: its contribution is subtracted explicitly.
       * Walk in original-X order, retaining original/native identity maps. */
      int lo = c->orderRank[i], hi = lo + 1;
      while (lo > 0) {
        const int id = c->order[lo - 1].id;
        const int native = int_TREE_X == NP_TREE_TRUE ? ipt_lookup_extern_X[id] : id;
        if (x->kw[native] == 0.0)
          break;
        --lo;
      }
      while (hi < n) {
        const int id = c->order[hi].id;
        const int native = int_TREE_X == NP_TREE_TRUE ? ipt_lookup_extern_X[id] : id;
        if (x->kw[native] == 0.0)
          break;
        ++hi;
      }
      /* One owning rank writes these exactly representable integer bounds.
       * Reuse the existing fold reduction; no additional collective. */
      a->folds[2 * n + i] = lo;
      a->folds[3 * n + i] = hi;
    }
    double denominator = 0.0, self = 0.0, sum = 0.0, error = 0.0;
    if(c->stable){
      NPLPSolveWorkspace *work=&x->regression_solve_workspace;
      NPConditionalLocalQRStatus rank=np_cqr_local_row(&a->qr.fast,work,n,1,terms,
        np_glp_cv_cache.terms,matrix_X_continuous_train_extern,x->kw,pos,x->mean_row,NULL);
      if(rank==NP_CQR_LOCAL_FULL){
        if(cqr_compensated_row(&a->qr.fast,n,terms,np_glp_cv_cache.terms,
          matrix_X_continuous_train_extern[0],x->kw,pos,c->coefficient_work,
          c->stable_coef+(size_t)i*terms,x->mean_row)){
          np_conditional_failure_record(NP_CONDITIONAL_NUMERICAL_FAILURE,i+1);
          fail=1;break;
        }
        /* Conversion must reproduce the independent QR row before its
         * coefficients enter compressed I1. Keep the frozen scaled row gate. */
        double row_scale=0.0,row_error=0.0;
        for(int j=0;j<n;j++) {
          const int original=int_TREE_X==NP_TREE_TRUE ? ipt_extern_X[j] : j;
          row_scale=fmax(row_scale,fabs(a->row[original]));
          row_error=fmax(row_error,fabs(x->mean_row[j]-a->row[original]));
          if(j!=pos && x->kw[j]>0)
            c->radius[i]=fmax(c->radius[i],fabs(matrix_X_continuous_train_extern[0][j]-matrix_X_continuous_train_extern[0][pos]));
        }
        if(!(row_error <= 1e-10*(1.0+row_scale))) {
          np_conditional_failure_record(NP_CONDITIONAL_COEFFICIENT_FAILURE,i+1);
          fail=1;break;
        }
      } else if(rank!=NP_CQR_LOCAL_DEFICIENT &&
                rank!=NP_CQR_LOCAL_AMBIGUOUS){
        np_conditional_failure_record(rank==NP_CQR_LOCAL_AMBIGUOUS ?
          NP_CONDITIONAL_RANK_AMBIGUOUS : NP_CONDITIONAL_NUMERICAL_FAILURE,i+1);
        fail=1;break;
      }
      denominator=1.;
    } else if (scalar) {
      for (int j = 0; j < n; ++j)
        denominator += x->kw[j];
    } else {
      for (int t = 0; t < terms; ++t)
        self += x->basis[t][pos] * x->regression_solve_workspace.rhs_work[t];
      self *= x->kw[pos];
      if (!np_lp_delete_denominator(self, &denominator)) {
        fail = 1;
        break;
      }
    }
    if (!R_FINITE(denominator) || denominator == 0.0) {
      fail = 1;
      break;
    }
    for (int t = 0; t < terms; ++t)
      c->coef[(size_t)terms*i+t] =
        /* FULL rows use stable_coef, not the unwritten incumbent solve RHS. */
        (c->stable && c->radius[i] > 0.0) ? 0.0 :
        (scalar ? 1.0 : x->regression_solve_workspace.rhs_work[t]) *
        (uniform && !c->stable ? k0 : 1.0) / denominator;
    /* Same one-owner transport: Epan-X now transports its canonical radius;
     * deletion is by disjoint ranges, not subtraction of the self term. */
    c->self[i] = c->stable ? h : uniform ? (scalar ? k0 : self) / denominator : h;
    for (int j = 0; j < n; ++j) {
      a->folds[n + i] += fabs(a->row[j]);
      np_cgnn_compensated(a->row[j] * a->yrow[j] / n, &sum, &error);
    }
    a->folds[i] = sum + error;
  }
#ifdef MPI2
  if (np_conditional_outer_buffer_finish(parallel, c->coefficient_count, fail, c->coef, NULL,
                                       "conditional GNN coefficients"))
    return R_NilValue;
  if (np_conditional_outer_buffer_finish(parallel, n, 0, c->self, NULL,
                                       "conditional GNN deletion coefficients"))
    return R_NilValue;
  if (np_conditional_outer_buffer_finish(parallel, foldPlanes * n, 0, a->folds, NULL,
                                       "conditional GNN I2 and influence bounds"))
    return R_NilValue;
#else
  if (fail)
    return R_NilValue;
#endif
  double cross = 0.0, cross_error = 0.0;
  for (int i = 0; i < n; ++i) {
    if (uniform) {
      c->first[i] = (int)a->folds[2 * n + i];
      c->end[i] = (int)a->folds[3 * n + i];
    } else {
      /* Raw distance/radius arithmetic matches the native support test and
       * does not collapse nearby observations in global normalized units. */
      const double h = c->self[i], xi = c->x[i];
      int left = 0, right = n;
      while (left < right) {
        const int mid = left+(right-left)/2;
        const double u = (c->order[mid].x-xi)/h;
        if (c->order[mid].x < xi && !(u*u < 5.0)) left = mid+1;
        else right = mid;
      }
      c->first[i] = left;
      right = n;
      while (left < right) {
        const int mid = left+(right-left)/2;
        const double u = (c->order[mid].x-xi)/h;
        if (u*u < 5.0) left = mid+1; else right = mid;
      }
      c->end[i] = left;
    }
    c->maxL1 = fmax(c->maxL1, a->folds[n + i]);
    np_cgnn_compensated(a->folds[i], &cross, &cross_error);
  }
  np_conditional_xrow_ctx_clear(x);
  np_conditional_yrow_ctx_clear(&a->yctx);
  if(c->stable){
    for(int i=0;i<n;i++)if(c->radius[i]>0)c->full_rows++;else c->deficient_rows++;
    fail=np_cgnn_local_prepare(&c->local,n,uniform?0:2,c->order,c->self,c->first,c->end);
    if(!fail && c->full_rows)fail=np_cgnn_stable_prepare(&c->stable_moments,n,
      terms-1+(uniform?0:2),c->order,c->self,c->first,c->end,c->radius);
    if(!fail && c->full_rows)fail=npc_gnn_recipes(&c->stable_moments,n,terms,c->order,
      c->radius,c->self,c->stable_coef,!uniform,k0);
  } else if (!uniform)
    fail = np_cgnn_local_prepare(&c->local, n, 2*(KERNEL_reg_extern-3),
                                  c->order, c->self, c->first, c->end);
  if(c->stable && fail)
    np_conditional_failure_record(NP_CONDITIONAL_NUMERICAL_FAILURE,0);
  first = 0;
  count = g->count;
#ifdef MPI2
  np_objective_outer_owned_rows(0, g->count, parallel, &first, &count);
#endif
  (void)np_mseries_accelerate_enabled();
#ifdef NP_CF167_TRACE
  REprintf("CF167 prefix rank=%d interval_first=%d interval_count=%d total=%d\n",
    np_cgnn_trace_rank(),first,count,g->count);
#endif
  int narrow_constants_ready = c->ykernel != 0;
  for (int at = first; at < first + count && !fail; ++at) {
    const NPGNNIntegralInterval interval = g->intervals[at];
    np_progress_bandwidth_loop_step();
    double sum = 0.0, error = 0.0;
    c->acceptedError = c->capped = 0.0;
    c->interval_visits=0;
    for (int slot = 0; slot < 2 && !fail; ++slot) {
      c->anchor = slot ? interval.successor_anchor : interval.primary_anchor;
      c->narrow = 0;
      int any = 0;
      for (int r = 0; r < n; ++r) {
        const int i = g->order[r];
        c->active[i] = ((r >= interval.first_deleted && r < interval.end_deleted) == slot);
        any += c->active[i];
      }
      if (!any)
        continue;
      c->activeCount = any;
      if (c->tree)
        np_cgnn_prepare_tree(c);
      /* The two zero-radius endpoint limits have opposite signs. */
      const double u0 = interval.lo == c->anchor ? R_PosInf : 1.0 / (interval.lo - c->anchor),
                   u1 = interval.hi == c->anchor ? R_NegInf : 1.0 / (interval.hi - c->anchor);
      double lo = fmin(u0, u1), hi = fmax(u0, u1);
      if (lo == hi && R_FINITE(interval.lo) && R_FINITE(interval.hi) &&
          interval.lo < interval.hi && interval.lo != c->anchor && interval.hi != c->anchor) {
        c->du = fabs((interval.hi - interval.lo) / (interval.lo - c->anchor) /
                     (interval.hi - c->anchor));
        if (!(c->du > 0.0) || !R_FINITE(c->du)) {
          fail = 1;
          break;
        }
        c->narrow = 1;
        c->qleft = interval.lo;
        lo = 0.0;
        hi = 1.0;
        if (!compact && !narrow_constants_ready) {
          np_cgnn_gaussian_derivative_constants(c->ykernel, c->logG);
          narrow_constants_ready = 1;
        }
      }
      if (!R_FINITE(lo) || !R_FINITE(hi)) {
        const int finite = np_cgnn_uniform_finite_panel(c, &lo, &hi);
        if (finite == 0)
          continue;
        if (finite < 0) {
          fail = 1;
          break;
        }
      }
      if (!(lo < hi)) {
        fail = 1;
        break;
      }
      const double before_error = c->acceptedError;
      if (compact)
        np_cgnn_compact_integrate(c, lo, hi, a->out);
      else
        np_cgnn_integrate(c, lo, hi, 0, a->out);
      if (c->failed) {
        fail = 1;
        break;
      }
      if (c->narrow) {
        for (int i = 0; i < n; ++i)
          a->out[i] *= c->du;
        c->acceptedError = before_error + (c->acceptedError - before_error) * c->du;
      }
      for (int i = 0; i < n; ++i)
        np_cgnn_compensated(a->out[i] / n, &sum, &error);
    }
    a->pieces[3 * at] = sum + error;
    a->pieces[3 * at + 1] = c->acceptedError;
    a->pieces[3 * at + 2] = c->capped;
  }
#ifdef MPI2
  if (np_conditional_outer_buffer_finish(parallel, 3 * g->count, fail, a->pieces, NULL,
                                       "conditional GNN response intervals"))
    return R_NilValue;
#else
  if (fail)
    return R_NilValue;
#endif
  double total = 0.0, correction = 0.0, estimate = 0.0, caps = 0.0;
  for (int at = 0; at < g->count; ++at) {
    np_cgnn_compensated(a->pieces[3 * at], &total, &correction);
    estimate += a->pieces[3 * at + 1];
    caps += a->pieces[3 * at + 2];
  }
  a->score = (total + correction) - 2.0 * (cross + cross_error);
#ifdef NP_CF167_TRACE
  REprintf("CF202 components I1=%.17g I2=%.17g visits=%zu full=%d deficient=%d\n",
    total+correction,cross+cross_error,c->visits,c->full_rows,c->deficient_rows);
#endif
  if (!R_FINITE(a->score) || !R_FINITE(estimate))
    return R_NilValue;
  np_cgnn_notice_caps += caps;
  if (caps > 0.0)
    np_cgnn_notice_error = fmax(np_cgnn_notice_error, estimate);
  a->status = 0;
  return R_NilValue;
}

static int np_cgnn_prefix_cvls(double *vsf, double *cv) {
  NPGNNConditionalCall call = {.vsf = vsf, .status = 1};
  R_UnwindProtect(np_cgnn_body, &call, np_cgnn_cleanup, &call, NULL);
  if (!call.status)
    *cv = call.score;
  return call.status;
}
