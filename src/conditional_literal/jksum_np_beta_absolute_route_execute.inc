static SEXP np_beta_absolute_route_execute(void *data)
{
  NPBetaAbsoluteRouteExecution * const execution =
    (NPBetaAbsoluteRouteExecution *)data;

  execution->status = np_beta_absolute_route_body(
    execution->call, &execution->owner);
  return R_NilValue;
}
