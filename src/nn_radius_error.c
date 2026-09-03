#include <stdio.h>
#include <R.h>
#include <Rinternals.h>
#include "headers.h"
#include "nn_radius_error.h"

static int np_nn_value_compare(const void *a, const void *b)
{
  const double x = *(const double *)a, y = *(const double *)b;
  return (x > y) - (x < y);
}

static int np_nn_equal_count(const double *sorted, int n, double value)
{
  int lo = 0, hi = n;
  while(lo < hi) {
    const int mid = lo + (hi - lo) / 2;
    if(sorted[mid] < value) lo = mid + 1;
    else hi = mid;
  }
  const int first = lo;
  hi = n;
  while(lo < hi) {
    const int mid = lo + (hi - lo) / 2;
    if(sorted[mid] <= value) lo = mid + 1;
    else hi = mid;
  }
  return lo - first;
}

/* Called only after a terminal owner has received ZERO_RADIUS. Counting
 * equal-valued donors explains that failure without sorting distances,
 * reconstructing weights, changing query geometry, or doing MPI work.
 * Sort a disposable column so explaining a failure is not quadratic in n.
 * Allocation failure retains the original cause, with less detail. */
NPNNZeroRadiusInfo np_nn_zero_radius_info(
  int bandwidth_type, int ntrain, int neval, int ncon,
  double **train, double **eval, double *scale,
  const NPNNGeometryContext *geometry, int coordinate_offset)
{
  NPNNZeroRadiusInfo info = {bandwidth_type, NA_INTEGER, NA_INTEGER,
    NA_INTEGER, NA_INTEGER, NA_REAL, NA_REAL};
  const int adaptive = bandwidth_type == BW_ADAP_NN;
  const NPNNQueryMode mode = adaptive ? NP_NN_QUERY_TRAINING_IDENTITY :
    (geometry == NULL ? NP_NN_QUERY_EXTERNAL : geometry->mode);

  /* Fold-specific support is not reconstructed by this explanation helper. */
  if(geometry != NULL && geometry->mode == NP_NN_QUERY_ADAPTIVE_FOLD_PREPARE)
    return info;
  if(ntrain < 1 || (size_t)ntrain > SIZE_MAX / sizeof(double))
    return info;
  double *sorted = (double *)malloc((size_t)ntrain * sizeof(double));
  if(sorted == NULL)
    return info;
  for(int d = 0; d < ncon; ++d) {
    int k, extended;
    double distance_scale;
    if(np_nn_lookup_from_scale(ntrain, 1, scale[d], &k,
                               &distance_scale, &extended) != 0)
      continue;
    int finite = 1;
    for(int j = 0; j < ntrain; ++j) {
      sorted[j] = train[d][j];
      if(!R_FINITE(sorted[j])) finite = 0;
    }
    if(!finite)
      continue;
    qsort(sorted, (size_t)ntrain, sizeof(double), np_nn_value_compare);
    for(int i = 0; i < (adaptive ? ntrain : neval); ++i) {
      const double value = adaptive ? train[d][i] : eval[d][i];
      if(!R_FINITE(value))
        continue;
      const int omit = mode == NP_NN_QUERY_TRAINING_IDENTITY ||
        (mode == NP_NN_QUERY_TRAINING_MAP && geometry != NULL &&
         geometry->eval_to_train != NULL && geometry->eval_to_train[i] >= 0);
      const int count = np_nn_equal_count(sorted, ntrain, value);
      if(count - omit >= k) {
        info.coordinate = coordinate_offset + d + 1;
        info.lookup_k = k;
        info.matching_donors = count - omit;
        info.excluded = omit;
        info.query_value = value;
        info.nn_index = scale[d];
        free(sorted);
        return info;
      }
    }
  }
  free(sorted);
  return info;
}

NPNNZeroRadiusInfo np_nn_zero_radius_info_flat(
  int bandwidth_type, int ntrain, int neval, int ncon,
  const double *train, const double *eval, const double *scale,
  const NPNNGeometryContext *geometry, int coordinate_offset)
{
  NPNNZeroRadiusInfo info = {bandwidth_type, NA_INTEGER, NA_INTEGER,
    NA_INTEGER, NA_INTEGER, NA_REAL, NA_REAL};
  for(int d = 0; d < ncon; ++d) {
    double *train_column = (double *)train + (size_t)d * ntrain;
    double *eval_column = eval == NULL ? train_column :
      (double *)eval + (size_t)d * neval;
    info = np_nn_zero_radius_info(bandwidth_type, ntrain, neval, 1,
      &train_column, &eval_column, (double *)scale + d, geometry, coordinate_offset + d);
    if(info.coordinate != NA_INTEGER)
      break;
  }
  return info;
}

void np_nn_zero_radius_error(const NPNNZeroRadiusInfo *info)
{
  const char *kind = info->bandwidth_type == BW_ADAP_NN ? "adaptive" : "generalized";
  char reason[512], message[576];
  if(info->coordinate != NA_INTEGER) {
    snprintf(reason, sizeof(reason),
      "repeated values give a zero literal radius for the %s nearest-neighbor bandwidth "
      "(NN index=%.17g, k=%d, query value=%.17g, %d matching donors; %s). "
      "Use fixed bandwidths (bwtype=\"fixed\", the default) instead.",
      kind, info->nn_index, info->lookup_k, info->query_value,
      info->matching_donors, info->excluded ? "one occurrence excluded" :
      "external query, no occurrence excluded");
    snprintf(message, sizeof(message), "continuous variable %d: %s",
             info->coordinate, reason);
  } else {
    snprintf(reason, sizeof(reason),
      "the %s nearest-neighbor bandwidth has a zero literal radius. "
      "Use fixed bandwidths (bwtype=\"fixed\", the default) instead.", kind);
    snprintf(message, sizeof(message), "%s", reason);
  }
  const char *fields[] = {"message", "call", "continuous.index", "bwtype",
    "nn.index", "lookup.k", "query.value", "matching.donors", "excluded", "reason"};
  SEXP condition = PROTECT(Rf_allocVector(VECSXP, 10));
  SEXP names = PROTECT(Rf_allocVector(STRSXP, 10));
  SEXP classes = PROTECT(Rf_allocVector(STRSXP, 4));
  for(int i = 0; i < 10; ++i)
    SET_STRING_ELT(names, i, Rf_mkChar(fields[i]));
  SET_VECTOR_ELT(condition, 0, Rf_mkString(message));
  SET_VECTOR_ELT(condition, 1, R_NilValue);
  SET_VECTOR_ELT(condition, 2, Rf_ScalarInteger(info->coordinate));
  SET_VECTOR_ELT(condition, 3, Rf_mkString(kind));
  SET_VECTOR_ELT(condition, 4, Rf_ScalarReal(info->nn_index));
  SET_VECTOR_ELT(condition, 5, Rf_ScalarInteger(info->lookup_k));
  SET_VECTOR_ELT(condition, 6, Rf_ScalarReal(info->query_value));
  SET_VECTOR_ELT(condition, 7, Rf_ScalarInteger(info->matching_donors));
  SET_VECTOR_ELT(condition, 8, Rf_ScalarInteger(info->excluded));
  SET_VECTOR_ELT(condition, 9, Rf_mkString(reason));
  Rf_setAttrib(condition, R_NamesSymbol, names);
  SET_STRING_ELT(classes, 0, Rf_mkChar("np_nn_zero_radius"));
  SET_STRING_ELT(classes, 1, Rf_mkChar("simpleError"));
  SET_STRING_ELT(classes, 2, Rf_mkChar("error"));
  SET_STRING_ELT(classes, 3, Rf_mkChar("condition"));
  Rf_classgets(condition, classes);
  /* Base stop preserves structured fields across existing native unwind
   * boundaries. Rf_error can transport only a string, not a typed condition. */
  SEXP call = PROTECT(Rf_lang2(Rf_install("stop"), condition));
  Rf_eval(call, R_BaseEnv);
  UNPROTECT(4);
  Rf_error("internal zero-radius condition unexpectedly returned");
}
