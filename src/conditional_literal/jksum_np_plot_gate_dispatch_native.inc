SEXP np_plot_gate_dispatch_native(SEXP kind, SEXP native_args, SEXP geometry,
  int support_columns, const int *support_counts,
  const double *const *support_values)
{
  NPPlotGateCall call = {0};
  if(TYPEOF(kind) != STRSXP || XLENGTH(kind) != 1 || STRING_ELT(kind, 0) == NA_STRING)
    error("plot gate dispatcher requires one fixed native operation");
  const char *label = CHAR(STRING_ELT(kind, 0));
  if(strcmp(label, "kernelsum") == 0) call.dispatch = 0;
  else if(strcmp(label, "lp_apply_ctx") == 0) call.dispatch = 1;
  else error("unknown plot gate native operation");
  call.args = native_args; call.geometry = geometry;
  call.support_columns = support_columns; call.parsed_counts = support_counts;
  call.parsed_support = support_values;
  return R_UnwindProtect(np_plot_gate_body, &call,
                         np_plot_gate_cleanup, &call, NULL);
}
