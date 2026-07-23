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

#endif
