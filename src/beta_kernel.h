#ifndef NP_BETA_KERNEL_H
#define NP_BETA_KERNEL_H

/*
 * Coordinate-aware associated beta-kernel primitives.
 *
 * These functions intentionally do not use the legacy unary continuous-kernel
 * registry.  A beta kernel needs the evaluation coordinate, observation
 * coordinate, finite support, and public metric bandwidth separately.
 */

typedef enum {
  NP_BETA_OK = 0,
  NP_BETA_ERR_NONFINITE = 1,
  NP_BETA_ERR_BOUNDS = 2,
  NP_BETA_ERR_BANDWIDTH = 3,
  NP_BETA_ERR_SCALE = 4,
  NP_BETA_ERR_OBSERVATION = 5,
  NP_BETA_ERR_CONCENTRATION = 6,
  NP_BETA_ERR_NUMERIC = 7,
  NP_BETA_ERR_RANGE = 8
} np_beta_status;

typedef enum {
  NP_BETA_EVAL_BELOW = -1,
  NP_BETA_EVAL_INSIDE = 0,
  NP_BETA_EVAL_ABOVE = 1
} np_beta_eval_location;

typedef struct {
  double lower;
  double upper;
  double support_length;
  double target_unit;
  double target_complement_unit;
  double observation_unit;
  double observation_complement_unit;
  double log_observation_unit;
  double log_observation_complement_unit;
  double concentration;
  np_beta_eval_location eval_location;
  int observation_endpoint;
} np_beta_shape;

np_beta_status np_beta_shape_init(double evaluation,
                                  double observation,
                                  double bandwidth,
                                  double lower,
                                  double upper,
                                  int scale,
                                  np_beta_shape *shape);

double np_beta_pdf_scale(const np_beta_shape *shape,
                         np_beta_status *status);

double np_beta_log_pdf_scale(const np_beta_shape *shape,
                             np_beta_status *status);

double np_beta_pdf_order2(double evaluation,
                          double observation,
                          double bandwidth,
                          double lower,
                          double upper,
                          np_beta_status *status);

double np_beta_pdf_order(double evaluation,
                         double observation,
                         double bandwidth,
                         double lower,
                         double upper,
                         int order,
                         np_beta_status *status);

double np_beta_log_pdf_order2(double evaluation,
                              double observation,
                              double bandwidth,
                              double lower,
                              double upper,
                              np_beta_status *status);

double np_beta_log_abs_pdf_order(double evaluation,
                                 double observation,
                                 double bandwidth,
                                 double lower,
                                 double upper,
                                 int order,
                                 int *sign,
                                 np_beta_status *status);

double np_beta_cdf_order2(double evaluation,
                          double observation,
                          double bandwidth,
                          double lower,
                          double upper,
                          np_beta_status *status);

double np_beta_cdf_order(double evaluation,
                         double observation,
                         double bandwidth,
                         double lower,
                         double upper,
                         int order,
                         np_beta_status *status);

double np_beta_log_overlap_order2(double center_one,
                                  double bandwidth_one,
                                  double center_two,
                                  double bandwidth_two,
                                  double lower,
                                  double upper,
                                  np_beta_status *status);

double np_beta_overlap_order2(double center_one,
                              double bandwidth_one,
                              double center_two,
                              double bandwidth_two,
                              double lower,
                              double upper,
                              np_beta_status *status);

double np_beta_overlap_order(double center_one,
                             double bandwidth_one,
                             double center_two,
                             double bandwidth_two,
                             double lower,
                             double upper,
                             int order,
                             np_beta_status *status);

int np_beta_order_supported(int order);

np_beta_status np_beta_signed_log_absolute(double positive_log,
                                           double negative_log,
                                           double *log_absolute,
                                           int *sign);

const char *np_beta_status_message(np_beta_status status);

#endif
