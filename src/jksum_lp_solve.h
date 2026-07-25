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
  int factor_ready;
  int factor_p;
} NPLPSolveWorkspace;

typedef struct {
  int n_capacity;
  int p_capacity;
  size_t xqr_capacity;
  double *xqr;
  double *qraux;
  double *work;
  double *y;
  double *qy;
  int *pivot;
} NPGLPQRDropWorkspace;

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
 * Exact basis-general influence row for a one-column weighted design:
 * w_i z_i z_eval / sum_j(w_j z_j^2).  positive_weights_only preserves the
 * QR owner's historical treatment of nonpositive weights; otherwise signed
 * higher-order kernel weights are retained.  output_stride permits both
 * contiguous native rows and strided R column-major matrix rows.  Returns
 * zero on success; it never allocates or calls BLAS/LAPACK.
 */
int np_lp_width_one_influence_row(const double *basis_train,
                                  int n,
                                  const double *kw,
                                  double basis_eval,
                                  double *row_out,
                                  size_t output_stride,
                                  int positive_weights_only);

/*
 * Reusable QR workspace for a leave-one-out local-polynomial influence row.
 * Width one computes the exact scalar influence row directly without QR
 * allocation or LAPACK.  Wider arithmetic and the dqrdc2/dqrqy transcript
 * match the historical helper; only allocation ownership moves from each row
 * to the enclosing owner.
 */
void np_glp_qr_drop_workspace_init(NPGLPQRDropWorkspace *workspace);
void np_glp_qr_drop_workspace_clear(NPGLPQRDropWorkspace *workspace);
int np_glp_qr_drop_workspace_reserve(NPGLPQRDropWorkspace *workspace,
                                     int n,
                                     int p);
int np_glp_qr_drop_workspace_apply(NPGLPQRDropWorkspace *workspace,
                                   double **basis,
                                   int n,
                                   int p,
                                   const double *kw,
                                   int eval_pos,
                                   double *row_out);

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
