static int np_conditional_lp_project_accepted(
  NPLPSolveWorkspace *workspace,
  const int nterms,
  const int nrhs,
  const double pristine_anchor,
  const NPLPSolvePolicyDiagnostics *diagnostics)
{
  const int ok = diagnostics->ridge_total > 0.0 ?
    np_lp_solve_workspace_solve_adjoint_factored(
      workspace, nterms, nrhs, pristine_anchor, diagnostics) == NP_LP_SOLVE_POLICY_OK :
    np_lp_solve_workspace_solve_factored(workspace, nterms, nrhs);
  if(ok && np_conditional_plot_capture_current != NULL)
    np_conditional_plot_capture_accepted(workspace, nterms, nrhs);
  return ok;
}
