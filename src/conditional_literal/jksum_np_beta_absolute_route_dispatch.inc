static int np_beta_absolute_route_dispatch(
  const NPBetaAbsoluteRouteCall * const call,
  const int suppress_parallel)
{
  if(call == NULL || (suppress_parallel != 0 && suppress_parallel != 1))
    return KWSNP_ERR_BADINVOC;
#ifdef MPI2
  if(iNum_Processors > 1 && !suppress_parallel &&
     !np_mpi_local_regression_active() &&
     call->regression_moment_context == NULL &&
     (call->bandwidth_mode == BW_FIXED ||
      call->bandwidth_mode == BW_GEN_NN))
    return np_beta_absolute_route_mpi_owner(call);
#endif
  return np_beta_absolute_route(call);
}
