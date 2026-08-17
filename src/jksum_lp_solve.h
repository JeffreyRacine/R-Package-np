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
  int *ipiv;
  double *rcond_work;
  int *rcond_iwork;
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
int np_lp_solve_workspace_solve(NPLPSolveWorkspace *workspace,
                                int p,
                                int nrhs);

/*
 * Canonical response-oriented bounded ridge transcript.  The caller fills
 * the pristine Gram and response-moment columns, and supplies the fixed ridge
 * increment owned by its statistical sample.  Failed Gram factorizations add that
 * increment to the source Gram diagonal in ascending order.  Once a ridged
 * system factors, the accumulated intercept restoration is applied to every
 * response RHS and all columns are solved from that retained factorization.
 * Response magnitude and batching therefore cannot select a different ridge.
 * The final factorization remains available to solve_factored().  This is
 * deliberately response-only: the
 * adjoint influence-row orientation has a different transform placement and
 * enters only with its own direct identity proof.
 */
NPLPSolvePolicyStatus np_lp_solve_workspace_solve_response(
  NPLPSolveWorkspace *workspace,
  int p,
  int nrhs,
  double ridge_increment,
  NPLPSolvePolicyDiagnostics *diagnostics);

/*
 * Adjoint sibling of solve_response().  It selects the identical Gram-owned
 * factor/ridge state, solves the evaluation-basis RHS columns with the
 * transpose of the retained LU, and only then applies the transposed
 * accumulated-ridge intercept transform to each solution.  The stored Gram is
 * mathematically symmetric, but independently accumulated moment columns need
 * not be bitwise symmetric; an ordinary non-transposed solve is therefore not
 * an adjoint.  Moving either transpose to the wrong side would not reproduce
 * the response fitted-value map.
 */
NPLPSolvePolicyStatus np_lp_solve_workspace_solve_adjoint(
  NPLPSolveWorkspace *workspace,
  int p,
  int nrhs,
  double ridge_increment,
  NPLPSolvePolicyDiagnostics *diagnostics);

/*
 * Solve adjoint RHS columns from the factor/ridge state retained by a prior
 * successful canonical response or adjoint solve.  The caller replaces only
 * rhs_source and supplies that retained solve's diagnostics.  This supports a
 * mixed response-plus-leverage consumer without refactorizing the same Gram or
 * pretending an evaluation-basis RHS is a response column.
 */
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
  size_t output_stride,
  NPLPSolvePolicyDiagnostics *diagnostics);

/*
 * Reusable contiguous Gram/RHS/rcond/solve storage for full-weight LP rows.
 * Width one uses scalar condition, solve, and inverse algebra and never
 * enters LAPACK.  The row owner reconstructs wider Gram and RHS buffers
 * before every call, so dgesv may overwrite them directly.  Wider systems
 * retain the same dsyev eigenvalue-ratio gate and dgesv transcript as the
 * historical row-fragmented route.
 */
void np_lp_full_row_workspace_init(NPLPFullRowWorkspace *workspace);
void np_lp_full_row_workspace_clear(NPLPFullRowWorkspace *workspace);
int np_lp_full_row_workspace_reserve(NPLPFullRowWorkspace *workspace,
                                     int p,
                                     int nrhs);
int np_lp_full_row_workspace_solve(NPLPFullRowWorkspace *workspace,
                                   int p,
                                   int nrhs,
                                   double min_rcond);

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
