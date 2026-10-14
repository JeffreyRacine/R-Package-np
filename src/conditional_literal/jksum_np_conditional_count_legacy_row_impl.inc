static NP_ALWAYS_INLINE NPContinuousKernelRowStatus
np_conditional_count_legacy_row_impl(const NPConditionalCountPlan * const plan,
                                const int is_x_side,
                                const int evaluation,
                                const NPBetaScaledRowCategoricalContext * const
                                  categorical_context,
                                double * const log_absolute,
                                signed char * const sign,
                                signed char * const categorical_sign,
                                double * const row,
                                double * const common_log_scale,
                                int * const bad_dimension,
                                const uint64_t * const support,
                                const uint64_t support_bit,
                                const double * const donor_mass,
                                double ** const donor_primary_bandwidth)
{
  const int dimensions = is_x_side ? plan->num_x_continuous :
    plan->num_y_continuous;
  const int categorical_dimensions = is_x_side ?
    plan->num_x_unordered + plan->num_x_ordered :
    plan->num_y_unordered + plan->num_y_ordered;
  const int do_cdf = is_x_side ? 0 : plan->do_distribution;
  const int bandwidth_mode = plan->bandwidth_mode;
  const int num_train = plan->num_train;
  double ** const train = is_x_side ? plan->train_x : plan->train_y;
  double ** const eval = is_x_side ? plan->eval_x : plan->eval_y;
  double ** const bandwidth_eval = is_x_side ?
    plan->bandwidth_eval_x : plan->bandwidth_eval_y;
  double ** const bandwidth_train = is_x_side ?
    plan->bandwidth_train_x : plan->bandwidth_train_y;
  const double * const lower = is_x_side ? plan->lower_x : plan->lower_y;
  const double * const upper = is_x_side ? plan->upper_x : plan->upper_y;
  const np_continuous_kernel_descriptor descriptor = is_x_side ?
    plan->descriptor_x : plan->descriptor_y;
  double maximum = -INFINITY;
  int observation;

  if(descriptor.family != NP_CKERNEL_FAMILY_LEGACY || dimensions < 0 ||
     categorical_dimensions < 0 ||
     dimensions > INT_MAX - categorical_dimensions ||
     dimensions + categorical_dimensions <= 0 ||
     (dimensions > 0 &&
      (train == NULL || eval == NULL || bandwidth_eval == NULL ||
       bandwidth_train == NULL || lower == NULL || upper == NULL)) ||
     log_absolute == NULL || sign == NULL || row == NULL ||
     common_log_scale == NULL ||
     (categorical_dimensions > 0 &&
      (categorical_context == NULL || categorical_sign == NULL)))
    return NP_CONTINUOUS_ROW_ERR_LAYOUT;

  for(observation = 0; observation < num_train; ++observation) {
    double product_log = 0.0;
    int product_sign = 1;
    int dimension;

    if((donor_mass != NULL && observation == evaluation) ||
       (support != NULL && !(support[observation] & support_bit))){
      log_absolute[observation]=-INFINITY;
      sign[observation]=0;
      continue;
    }

    for(dimension = 0; dimension < dimensions; ++dimension) {
      const double bandwidth = bandwidth_mode == BW_FIXED ?
        bandwidth_eval[dimension][0] :
        (bandwidth_mode == BW_GEN_NN ?
          bandwidth_eval[dimension][evaluation] :
          bandwidth_train[dimension][observation]);
      double scalar_log = -INFINITY;
      int scalar_sign = 0;
      const NPContinuousKernelScalarStatus scalar_status =
        donor_mass != NULL ?
        np_continuous_kernel_legacy_pdf_log_mass(
          descriptor.legacy_code, train[dimension][observation],
          eval[dimension][evaluation], bandwidth,
          donor_mass[2U*((size_t)dimension*(size_t)num_train +
                         (size_t)observation) +
                     (bandwidth != donor_primary_bandwidth[dimension][observation])],
          &scalar_log, &scalar_sign) :
        np_continuous_kernel_scalar_log(
          descriptor.family, descriptor.legacy_code, descriptor.order,
          do_cdf, eval[dimension][evaluation], train[dimension][observation],
          bandwidth, lower[dimension], upper[dimension],
          &scalar_log, &scalar_sign);

      if(scalar_status == NP_CONTINUOUS_KERNEL_SCALAR_ERR_KERNEL) {
        if(bad_dimension != NULL)
          *bad_dimension = (is_x_side ? 0 : plan->num_x_continuous) + dimension;
        return NP_CONTINUOUS_ROW_ERR_KERNEL;
      }
      if(scalar_status != NP_CONTINUOUS_KERNEL_SCALAR_OK) {
        if(bad_dimension != NULL)
          *bad_dimension = (is_x_side ? 0 : plan->num_x_continuous) + dimension;
        return scalar_status == NP_CONTINUOUS_KERNEL_SCALAR_ERR_NUMERIC ?
          NP_CONTINUOUS_ROW_ERR_NUMERIC : NP_CONTINUOUS_ROW_ERR_LAYOUT;
      }
      if(scalar_sign == 0) {
        product_log = -INFINITY;
        product_sign = 0;
        break;
      }
      product_log += scalar_log;
      product_sign *= scalar_sign;
      if(ISNAN(product_log) || product_log == INFINITY) {
        if(bad_dimension != NULL)
          *bad_dimension = (is_x_side ? 0 : plan->num_x_continuous) + dimension;
        return NP_CONTINUOUS_ROW_ERR_NUMERIC;
      }
    }
    log_absolute[observation] = product_log;
    sign[observation] = (signed char)product_sign;
    if(product_sign != 0 && product_log > maximum)
      maximum = product_log;
  }

  if(categorical_dimensions > 0) {
    const NPContinuousKernelRowStatus categorical_status =
      np_beta_categorical_log_factor(
        categorical_context, evaluation, -1, num_train,
        row, categorical_sign);

    if(categorical_status != NP_CONTINUOUS_ROW_OK)
      return categorical_status;
    maximum = -INFINITY;
    for(observation = 0; observation < num_train; ++observation) {
      if(sign[observation] == 0 || categorical_sign[observation] == 0) {
        log_absolute[observation] = -INFINITY;
        sign[observation] = 0;
      } else {
        log_absolute[observation] += row[observation];
        sign[observation] = (signed char)(
          (int)sign[observation] * (int)categorical_sign[observation]);
        if(ISNAN(log_absolute[observation]) ||
           log_absolute[observation] == INFINITY)
          return NP_CONTINUOUS_ROW_ERR_NUMERIC;
        if(log_absolute[observation] > maximum)
          maximum = log_absolute[observation];
      }
    }
  }

  if(maximum == -INFINITY) {
    for(observation = 0; observation < num_train; ++observation)
      row[observation] = 0.0;
    *common_log_scale = -INFINITY;
    return NP_CONTINUOUS_ROW_OK;
  }

  for(observation = 0; observation < num_train; ++observation) {
    row[observation] = sign[observation] == 0 ? 0.0 :
      (double)sign[observation] * exp(log_absolute[observation] - maximum);
    if(!R_FINITE(row[observation]))
      return NP_CONTINUOUS_ROW_ERR_NUMERIC;
  }
  *common_log_scale = maximum;
  return NP_CONTINUOUS_ROW_OK;
}
