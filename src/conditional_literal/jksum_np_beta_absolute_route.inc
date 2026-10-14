static int np_beta_absolute_route(const NPBetaAbsoluteRouteCall *call)
{
  NPBetaAbsoluteRouteExecution execution;

  if(call == NULL)
    return KWSNP_ERR_BADINVOC;
  execution.call = call;
  execution.status = KWSNP_ERR_BADINVOC;
  np_beta_absolute_route_owner_init(&execution.owner);
  execution.owner.route_row = call->row;
  execution.owner.owned_bandwidth_tmatrix =
    call->owned_bandwidth_tmatrix;
  R_UnwindProtect(
    np_beta_absolute_route_execute, &execution,
    np_beta_absolute_route_owner_cleanup, &execution.owner, NULL);
  return execution.status;
}
