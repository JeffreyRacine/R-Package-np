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

/*
 * The caller owns the workspace and its lifetime.  gram_source/rhs_source are
 * the caller-mutable pristine system; solve copies them to gram_work/rhs_work,
 * lets LAPACK overwrite only the work buffers, and returns the solution in
 * rhs_work.  This preserves the source system for ridge retries and additional
 * right-hand sides without process-global scratch or per-solve allocation.
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
 * Reusable QR workspace for a leave-one-out local-polynomial influence row.
 * The arithmetic and dqrdc2/dqrqy transcript match the historical helper;
 * only allocation ownership moves from each row to the enclosing owner.
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
 * The row owner reconstructs Gram and RHS before every call, so dgesv may
 * overwrite them directly.  The rcond gate uses the same dsyev eigenvalue
 * ratio and the solve uses the same dgesv transcript as the historical
 * MATRIX-based route.
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
 * to the preserved source.  This matches legacy mat_inv() retry ownership
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
