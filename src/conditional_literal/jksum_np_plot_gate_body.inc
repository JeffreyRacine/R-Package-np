static SEXP np_plot_gate_body(void *data)
{
  NPPlotGateCall *c = (NPPlotGateCall *)data;
  np_plot_gate_prepare(c);
  c->saved = np_gate_override_ctx;
  c->installed = 1;
  np_gate_override_ctx = c->gate;
  /* Fixed native dispatch only: no eval(), R closure, or callback argument. */
#define A(i) VECTOR_ELT(c->args, (i))
  if(c->dispatch == 0)
    return C_np_kernelsum(A(0), A(1), A(2), A(3), A(4), A(5), A(6), A(7),
      A(8), A(9), A(10), A(11), A(12), A(13), A(14), A(15), A(16), A(17), A(18));
  return C_np_regression_lp_apply_conditional_ctx(A(0), A(1), A(2), A(3),
    A(4), A(5), A(6), A(7), A(8), A(9), A(10), A(11), A(12), A(13),
    A(14), A(15), A(16), A(17), A(18), A(19), A(20), A(21), A(22),
    A(23), A(24), A(25), A(26));
#undef A
}
