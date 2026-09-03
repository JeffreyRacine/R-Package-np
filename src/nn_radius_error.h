#ifndef NP_NN_RADIUS_ERROR_H
#define NP_NN_RADIUS_ERROR_H

/* Value-only, failure-path diagnostic; never retained across calls/ranks. */
typedef struct {
  int bandwidth_type;
  int coordinate;
  int lookup_k;
  int matching_donors;
  int excluded;
  double query_value;
  double nn_index;
} NPNNZeroRadiusInfo;

NPNNZeroRadiusInfo np_nn_zero_radius_info(
  int bandwidth_type, int ntrain, int neval, int ncon,
  double **train, double **eval, double *scale,
  const NPNNGeometryContext *geometry, int coordinate_offset);
void np_nn_zero_radius_error(const NPNNZeroRadiusInfo *info);
NPNNZeroRadiusInfo np_nn_zero_radius_info_flat(
  int bandwidth_type, int ntrain, int neval, int ncon,
  const double *train, const double *eval, const double *scale,
  const NPNNGeometryContext *geometry, int coordinate_offset);

#endif
