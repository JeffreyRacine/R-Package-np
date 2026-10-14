static NPContinuousKernelRowStatus
np_conditional_count_legacy_row(const NPConditionalCountPlan * const plan,
                                const int is_x_side,
                                const int evaluation,
                                const NPBetaScaledRowCategoricalContext * const
                                  categorical_context,
                                double * const log_absolute,
                                signed char * const sign,
                                signed char * const categorical_sign,
                                double * const row,
                                double * const common_log_scale,
                                int * const bad_dimension)
{
  return np_conditional_count_legacy_row_impl(
    plan,is_x_side,evaluation,categorical_context,log_absolute,sign,
    categorical_sign,row,common_log_scale,bad_dimension,NULL,0,NULL,NULL);
}
