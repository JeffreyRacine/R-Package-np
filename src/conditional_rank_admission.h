#ifndef NP_CONDITIONAL_RANK_ADMISSION_H
#define NP_CONDITIONAL_RANK_ADMISSION_H
/* Cold classifier arithmetic copied unchanged from jksum_lp_solve.c. */
static int np_conditional_rank_scratch_shape(
  const NPLPSolveWorkspace *workspace, int p)
{
  return (workspace != NULL) && (p > 0) &&
    (workspace->p_capacity >= p) &&
    (workspace->gram_source != NULL) &&
    (workspace->gram_work != NULL) &&
    (workspace->rank_values != NULL) &&
    (workspace->rank_work != NULL) &&
    (workspace->rank_work_capacity >= 5U*(size_t)p);
}

static NPLPFactorAdmissionStatus np_conditional_cold_rank(
  NPLPSolveWorkspace *workspace, int p)
{
  const char jobu = 'N';
  const char jobvt = 'N';
  const int ldu = 1;
  const int ldvt = 1;
  const int lwork = 5*p;
  const double pu = (double)p*DBL_EPSILON;
  const double gamma_p = (pu < 1.0) ? pu/(1.0 - pu) : 1.0;
  double dummy_u = 0.0;
  double dummy_vt = 0.0;
  double max_singular;
  double min_singular;
  int info = 0;
  int i, j;

  if(!np_conditional_rank_scratch_shape(workspace, p))
    return NP_LP_FACTOR_ADMISSION_INVALID;

  workspace->factor_ready = 0;
  workspace->factor_p = 0;
  for(i = 0; i < p; i++){
    const double diagonal = fabs(workspace->gram_source[i + i*p]);
    double reference = diagonal;

    if(!R_FINITE(diagonal))
      return NP_LP_FACTOR_ADMISSION_NONFINITE;
    if(!(reference > 0.0)){
      for(j = 0; j < p; j++){
        const double value = fabs(workspace->gram_source[i + j*p]);
        if(!R_FINITE(value))
          return NP_LP_FACTOR_ADMISSION_NONFINITE;
        if(value > reference)
          reference = value;
      }
    }
    if(!(reference > 0.0))
      return NP_LP_FACTOR_ADMISSION_RANK_DEFICIENT;
    workspace->rank_values[i] = 1.0/sqrt(reference);
    if(!R_FINITE(workspace->rank_values[i]))
      return NP_LP_FACTOR_ADMISSION_NONFINITE;
  }

  for(j = 0; j < p; j++){
    for(i = 0; i < p; i++){
      workspace->gram_work[i + j*p] = workspace->gram_source[i + j*p]*
        workspace->rank_values[i]*workspace->rank_values[j];
      if(!R_FINITE(workspace->gram_work[i + j*p]))
        return NP_LP_FACTOR_ADMISSION_NONFINITE;
    }
  }

  F77_CALL(dgesvd)(&jobu, &jobvt, &p, &p,
                   workspace->gram_work, &p,
                   workspace->rank_values,
                   &dummy_u, &ldu,
                   &dummy_vt, &ldvt,
                   workspace->rank_work, &lwork,
                   &info FCONE FCONE);
  if(info != 0)
    return NP_LP_FACTOR_ADMISSION_FAILED;
  max_singular = workspace->rank_values[0];
  min_singular = workspace->rank_values[p - 1];
  if(!R_FINITE(max_singular) || !R_FINITE(min_singular))
    return NP_LP_FACTOR_ADMISSION_NONFINITE;
  if(!(max_singular > 0.0) ||
     (min_singular <= gamma_p*max_singular))
    return NP_LP_FACTOR_ADMISSION_RANK_DEFICIENT;
  return NP_LP_FACTOR_ADMISSION_REFACTOR;
}
/* Conditional positive-weight deleted designs must pass the authoritative
 * spectral rank rule. LU pivots are not a sufficient full-rank certificate.
 * Keep the shared regression solver and its ridge amount/correction intact. */
static NPLPSolvePolicyStatus np_conditional_solve_adjoint_ranked(
  NPLPSolveWorkspace *workspace, int p, int nrhs, double ridge_fraction,
  int rank_upper_bound, NPLPSolvePolicyDiagnostics *diagnostics)
{
  if(rank_upper_bound < NP_LP_RANK_UPPER_BOUND_UNKNOWN)
    return NP_LP_SOLVE_POLICY_INVALID;
  if(p > 1 && (rank_upper_bound == NP_LP_RANK_UPPER_BOUND_UNKNOWN ||
               rank_upper_bound >= p)){
    const NPLPFactorAdmissionStatus admission =
      np_conditional_cold_rank(workspace, p);
    if(admission == NP_LP_FACTOR_ADMISSION_RANK_DEFICIENT)
      rank_upper_bound = 0;
    else if(admission == NP_LP_FACTOR_ADMISSION_NONFINITE)
      return NP_LP_SOLVE_POLICY_NONFINITE;
    else if(admission == NP_LP_FACTOR_ADMISSION_INVALID)
      return NP_LP_SOLVE_POLICY_INVALID;
    else if(admission != NP_LP_FACTOR_ADMISSION_REFACTOR)
      return NP_LP_SOLVE_POLICY_FINAL_FAILED;
  }
  return np_lp_solve_workspace_solve_adjoint_ranked(workspace, p, nrhs,
    ridge_fraction, rank_upper_bound, diagnostics);
}

#endif
