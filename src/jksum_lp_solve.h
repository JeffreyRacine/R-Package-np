#ifndef NP_JKSUM_LP_SOLVE_H
#define NP_JKSUM_LP_SOLVE_H

#include <stddef.h>

typedef struct {
  int p_capacity;
  int nrhs_capacity;
  size_t gram_capacity;
  size_t rhs_capacity;
  double *gram_source;
  double *rhs_source;
  double *gram_work;
  double *rhs_work;
  double *rank_values;
  double *rank_work;
  size_t rank_work_capacity;
  int *ipiv;
  int factor_ready;
  int factor_p;
} NPLPSolveWorkspace;

typedef enum {
  NP_LP_SOLVE_POLICY_OK = 0,
  NP_LP_SOLVE_POLICY_INVALID,
  NP_LP_SOLVE_POLICY_NONFINITE,
  NP_LP_SOLVE_POLICY_RIDGE_EXHAUSTED,
  NP_LP_SOLVE_POLICY_FINAL_FAILED
} NPLPSolvePolicyStatus;

typedef struct {
  int ridge_steps;
  double ridge_total;
} NPLPSolvePolicyDiagnostics;

/*
 * A successful LU can still contain a working-precision zero pivot.  The
 * ordinary case retains its factor; only the O(p) diagonal screen's narrow
 * ambiguity band enters the cold, equilibrated singular-value adjudicator.
 * A cold full-rank result has overwritten the LU and therefore explicitly
 * requests one refactorization from its caller.
 */
typedef enum {
  NP_LP_FACTOR_ADMISSION_RETAINED = 0,
  NP_LP_FACTOR_ADMISSION_REFACTOR,
  NP_LP_FACTOR_ADMISSION_RANK_DEFICIENT,
  NP_LP_FACTOR_ADMISSION_NONFINITE,
  NP_LP_FACTOR_ADMISSION_INVALID,
  NP_LP_FACTOR_ADMISSION_FAILED
} NPLPFactorAdmissionStatus;

#define NP_LP_RANK_UPPER_BOUND_UNKNOWN (-1)

typedef struct {
  int p_capacity;
  int nrhs_capacity;
  size_t gram_capacity;
  size_t rhs_capacity;
  int rcond_lwork_capacity;
  int inverse_lwork_capacity;
  double *gram;
  double *rhs;
  int *ipiv;
  double *matrix_copy;
  double *rcond_values;
  double *rcond_work;
  double *inverse_work;
} NPLPFullRowWorkspace;

#define NP_LP_SOLVE_MAX_RIDGE_STEPS 128

typedef enum {
  NP_LP_WIDTH_ONE_OK = 0,
  NP_LP_WIDTH_ONE_INVALID,
  NP_LP_WIDTH_ONE_NONFINITE,
  NP_LP_WIDTH_ONE_RIDGE_FAILED
} NPLPWidthOneStatus;

/*
 * The caller owns the workspace and its lifetime.  gram_source/rhs_source are
 * the caller-mutable pristine system.  Width one is solved directly as scalar
 * algebra; it never enters LAPACK.  Wider systems copy into
 * gram_work/rhs_work and let LAPACK overwrite only the work buffers.  Both
 * paths return the solution in rhs_work and preserve the source system for
 * ridge retries and additional right-hand sides without process-global
 * scratch or per-solve allocation.
 */
void np_lp_solve_workspace_init(NPLPSolveWorkspace *workspace);
void np_lp_solve_workspace_clear(NPLPSolveWorkspace *workspace);
int np_lp_solve_workspace_reserve(NPLPSolveWorkspace *workspace,
                                  int p,
                                  int nrhs);

/* Reuse a factor/ridge state retained by a canonical response/adjoint solve. */
NPLPSolvePolicyStatus np_lp_solve_workspace_solve_adjoint_factored(
  NPLPSolveWorkspace *workspace,
  int p,
  int nrhs,
  const NPLPSolvePolicyDiagnostics *diagnostics);
/* Cold-path validation used only after an ordinary solve has failed. */
#if defined(__GNUC__) || defined(__clang__)
#define NP_LP_COLD_PATH __attribute__((cold))
#else
#define NP_LP_COLD_PATH
#endif
int np_lp_solve_workspace_sources_finite(
  const NPLPSolveWorkspace *workspace,
  int p,
  int nrhs) NP_LP_COLD_PATH;
#undef NP_LP_COLD_PATH

/*
 * Solve new pristine RHS columns with the LU/pivots retained by the most
 * recent successful solve.  This never refactorizes or changes the source
 * Gram and is valid only until reserve, clear, or another solve.
 */
int np_lp_solve_workspace_solve_factored(NPLPSolveWorkspace *workspace,
                                         int p,
                                         int nrhs);

/*
 * Canonical rank/admission primitives.  They are internal C interfaces, not
 * registered R entry points.  rank_upper_bound may be UNKNOWN; otherwise it
 * is an exact structural upper bound supplied by an accumulation owner.
 */
NPLPFactorAdmissionStatus np_lp_solve_workspace_admit_factor(
  NPLPSolveWorkspace *workspace,
  int p,
  int rank_upper_bound);
int np_lp_solve_workspace_ridge_increment(
  const NPLPSolveWorkspace *workspace,
  int p,
  double ridge_fraction,
  double *ridge_increment_out);

/*
 * Canonical policy entry points.  ridge_fraction is converted once from the
 * pristine Gram into an owner-invariant ridge increment.  rank_upper_bound
 * may be UNKNOWN; otherwise it is the exact number of nonzero donor rows and
 * therefore a structural upper bound on rank before regularization. Ordinary
 * response rows retain one-call DGESV and the successful LU for later
 * RHS/adjoint reuse.
 */
NPLPSolvePolicyStatus np_lp_solve_workspace_solve_response_ranked(
  NPLPSolveWorkspace *workspace,
  int p,
  int nrhs,
  double ridge_fraction,
  int rank_upper_bound,
  NPLPSolvePolicyDiagnostics *diagnostics);
NPLPSolvePolicyStatus np_lp_solve_workspace_solve_adjoint_ranked(
  NPLPSolveWorkspace *workspace,
  int p,
  int nrhs,
  double ridge_fraction,
  int rank_upper_bound,
  NPLPSolvePolicyDiagnostics *diagnostics);

/*
 * Exact structural rank upper bound for one weighted design row.  Each finite
 * or non-finite nonzero weight can contribute at most one independent donor
 * row, so the count is capped at p and requires no allocation.  Numerical
 * finiteness remains the solve policy's responsibility after accumulation.
 */
int np_lp_rank_upper_bound_from_weights(const double *weights, int n, int p);

/*
 * Exact basis-general influence row for a one-column signed weighted design:
 * w_i z_i z_eval / sum_j(w_j z_j^2).  output_stride permits both contiguous
 * native rows and strided R column-major matrix rows.  A finite zero
 * denominator follows the canonical bounded scalar ridge transcript used by
 * the hat owner; a nonzero finite denominator is never thresholded or
 * perturbed.  Returns NP_LP_WIDTH_ONE_OK on success and a typed failure
 * status otherwise.  On failure, row_out contents are undefined and must not
 * be consumed.  It never allocates or calls BLAS/LAPACK.
 */
NPLPWidthOneStatus np_lp_width_one_influence_row(
  const double *basis_train,
  int n,
  const double *kw,
  double basis_eval,
  double *row_out,
  size_t output_stride);

/*
 * Reusable contiguous storage for retained full-Gram inverse owners. Width
 * one uses scalar condition and inverse algebra and never enters LAPACK.
 * Direct response and influence solves use NPLPSolveWorkspace above; this
 * older workspace remains only for the all-row inverse topology.
 */
void np_lp_full_row_workspace_init(NPLPFullRowWorkspace *workspace);
void np_lp_full_row_workspace_clear(NPLPFullRowWorkspace *workspace);
int np_lp_full_row_workspace_reserve(NPLPFullRowWorkspace *workspace,
                                     int p,
                                     int nrhs);

/*
 * Invert the symmetric Gram buffer in place after the same dsyev rcond gate.
 * The retained inverse remains in gram; pivot, eigen, and dgetri work storage
 * are workspace-owned and reusable.
 */
int np_lp_full_row_workspace_invert(NPLPFullRowWorkspace *workspace,
                                    int p,
                                    double min_rcond);

/*
 * Retry an ungated retained inverse while preserving the source Gram in
 * matrix_copy.  The caller fills matrix_copy; each attempt copies it into the
 * destructive LAPACK buffer, and a failed attempt adds one fixed ridge step
 * to the preserved source. This preserves the historical retry ownership
 * without row-fragmented matrices or per-attempt allocation.
 */
int np_lp_full_row_workspace_invert_retryable(
  NPLPFullRowWorkspace *workspace,
  int p,
  double ridge_increment,
  int max_ridge_steps);

/*
 * After a successful retained inversion, reuse matrix_copy as a row-major
 * inverse view for row-oriented fitted-value and quadratic-form consumers.
 * No additional p-by-p storage is allocated.
 */
int np_lp_full_row_workspace_pack_inverse_rows(
  NPLPFullRowWorkspace *workspace,
  int p);

#endif
