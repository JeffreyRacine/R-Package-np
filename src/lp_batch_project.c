#include <limits.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>

#include <R.h>
#include <Rinternals.h>

#include "jksum_lp_solve.h"

/*
 * Internal fixed-LP bootstrap bridge.
 *
 * packed_gram is nsystem x p(p+1)/2 in upper-triangle order
 * (0,0), (0,1), ..., (p-1,p-1).  moments is either an nsystem x p
 * matrix for one response or a length-p list of nsystem x nrhs matrices.
 * projection has length p and represented_mass has one positive finite entry
 * per system.  The bridge owns only marshaling and projection; all condition,
 * ridge, factorization, and response-solve policy remains in
 * np_lp_solve_workspace_solve_response().
 *
 * The return value is a five-element list: values, status, failed_system,
 * ridge_steps, ridge_total.  values is committed only after the complete batch
 * succeeds.  The diagnostic vectors are NULL unless requested.  status uses
 * NPLPSolvePolicyStatus and failed_system is one-based (zero on success).
 */

static int np_lp_batch_matrix_dims(SEXP x, int *nrow, int *ncol)
{
  SEXP dim = getAttrib(x, R_DimSymbol);
  R_xlen_t expected = 0;

  if((TYPEOF(dim) != INTSXP) || (XLENGTH(dim) != 2))
    return 0;
  *nrow = INTEGER(dim)[0];
  *ncol = INTEGER(dim)[1];
  if((*nrow < 0) || (*ncol < 0) ||
     ((*ncol > 0) && ((R_xlen_t)*nrow >
                      R_XLEN_T_MAX/(R_xlen_t)*ncol)))
    return 0;
  expected = (R_xlen_t)*nrow*(R_xlen_t)*ncol;
  return XLENGTH(x) == expected;
}

static SEXP np_lp_batch_result(SEXP values,
                               int status,
                               int failed_system,
                               SEXP ridge_steps,
                               SEXP ridge_total)
{
  SEXP out = PROTECT(allocVector(VECSXP, 5));
  SEXP names = PROTECT(allocVector(STRSXP, 5));

  SET_STRING_ELT(names, 0, mkChar("values"));
  SET_STRING_ELT(names, 1, mkChar("status"));
  SET_STRING_ELT(names, 2, mkChar("failed_system"));
  SET_STRING_ELT(names, 3, mkChar("ridge_steps"));
  SET_STRING_ELT(names, 4, mkChar("ridge_total"));
  SET_VECTOR_ELT(out, 0, values);
  SET_VECTOR_ELT(out, 1, ScalarInteger(status));
  SET_VECTOR_ELT(out, 2, ScalarInteger(failed_system));
  SET_VECTOR_ELT(out, 3, ridge_steps);
  SET_VECTOR_ELT(out, 4, ridge_total);
  setAttrib(out, R_NamesSymbol, names);

  UNPROTECT(2);
  return out;
}

SEXP C_np_lp_batch_project(SEXP packed_gram,
                           SEXP moments,
                           SEXP projection,
                           SEXP represented_mass,
                           SEXP diagnostics)
{
  int nsystem = 0;
  int packed_width = 0;
  int p = 0;
  int nrhs = 0;
  int moment_nrow = 0;
  int moment_ncol = 0;
  int want_diagnostics = 0;
  int status = NP_LP_SOLVE_POLICY_OK;
  int failed_system = 0;
  int nprotect = 0;
  size_t expected_packed = 0;
  const double **moment_ptr = NULL;
  const double *packed_ptr = NULL;
  const double *projection_ptr = NULL;
  const double *mass_ptr = NULL;
  double *value_ptr = NULL;
  NPLPSolveWorkspace workspace;
  SEXP values = R_NilValue;
  SEXP ridge_steps = R_NilValue;
  SEXP ridge_total = R_NilValue;
  SEXP result = R_NilValue;

  if((TYPEOF(packed_gram) != REALSXP) ||
     !np_lp_batch_matrix_dims(packed_gram, &nsystem, &packed_width) ||
     (nsystem <= 0))
    error("packed_gram must be a non-empty double matrix");
  if((TYPEOF(projection) != REALSXP) ||
     (XLENGTH(projection) <= 0) || (XLENGTH(projection) > INT_MAX))
    error("projection must be a non-empty double vector");
  p = (int)XLENGTH(projection);
  if(((size_t)p > SIZE_MAX/((size_t)p + 1U)) ||
     ((size_t)p*((size_t)p + 1U)/2U > (size_t)INT_MAX))
    error("projection width is too large");
  expected_packed = (size_t)p*((size_t)p + 1U)/2U;
  if((size_t)packed_width != expected_packed)
    error("packed_gram width does not match projection length");
  if((TYPEOF(represented_mass) != REALSXP) ||
     (XLENGTH(represented_mass) != (R_xlen_t)nsystem))
    error("represented_mass must be a double vector with one entry per system");
  if((TYPEOF(diagnostics) != LGLSXP) || (XLENGTH(diagnostics) != 1) ||
     (LOGICAL(diagnostics)[0] == NA_LOGICAL))
    error("diagnostics must be TRUE or FALSE");
  want_diagnostics = LOGICAL(diagnostics)[0];

  if(TYPEOF(moments) == REALSXP) {
    if(!np_lp_batch_matrix_dims(moments, &moment_nrow, &moment_ncol) ||
       (moment_nrow != nsystem) || (moment_ncol != p))
      error("single-response moments must be an nsystem by p double matrix");
    nrhs = 1;
  } else if(TYPEOF(moments) == VECSXP) {
    if(XLENGTH(moments) != (R_xlen_t)p)
      error("multi-response moments must be a length-p list");
    for(int term = 0; term < p; term++) {
      SEXP component = VECTOR_ELT(moments, term);
      int component_nrow = 0;
      int component_ncol = 0;

      if((TYPEOF(component) != REALSXP) ||
         !np_lp_batch_matrix_dims(component,
                                  &component_nrow,
                                  &component_ncol) ||
         (component_nrow != nsystem) || (component_ncol <= 0) ||
         ((term > 0) && (component_ncol != nrhs)))
        error("multi-response moment components must be conformable double matrices");
      if(term == 0)
        nrhs = component_ncol;
    }
  } else {
    error("moments must be a double matrix or a list of double matrices");
  }

  if((nrhs <= 0) || ((size_t)nsystem > SIZE_MAX/(size_t)nrhs) ||
     ((R_xlen_t)nsystem > R_XLEN_T_MAX/(R_xlen_t)nrhs))
    error("batch result dimensions are too large");

  for(int term = 0; term < p; term++)
    if(!R_FINITE(REAL(projection)[term]))
      error("projection contains a non-finite value");
  for(int system = 0; system < nsystem; system++)
    if(!R_FINITE(REAL(represented_mass)[system]) ||
       (REAL(represented_mass)[system] <= 0.0))
      error("represented_mass must contain positive finite values");

  values = PROTECT(allocMatrix(REALSXP, nsystem, nrhs));
  nprotect++;
  value_ptr = REAL(values);
  if(want_diagnostics) {
    ridge_steps = PROTECT(allocVector(INTSXP, nsystem));
    nprotect++;
    ridge_total = PROTECT(allocVector(REALSXP, nsystem));
    nprotect++;
    for(int system = 0; system < nsystem; system++) {
      INTEGER(ridge_steps)[system] = NA_INTEGER;
      REAL(ridge_total)[system] = NA_REAL;
    }
  }

  if(TYPEOF(moments) == VECSXP) {
    if((size_t)p > SIZE_MAX/sizeof(*moment_ptr)) {
      UNPROTECT(nprotect);
      error("multi-response moment pointer table is too large");
    }
    moment_ptr = (const double **)R_alloc((size_t)p, sizeof(*moment_ptr));
    for(int term = 0; term < p; term++)
      moment_ptr[term] = REAL(VECTOR_ELT(moments, term));
  }

  packed_ptr = REAL(packed_gram);
  projection_ptr = REAL(projection);
  mass_ptr = REAL(represented_mass);
  np_lp_solve_workspace_init(&workspace);
  if(!np_lp_solve_workspace_reserve(&workspace, p, nrhs)) {
    np_lp_solve_workspace_clear(&workspace);
    UNPROTECT(nprotect);
    error("unable to allocate canonical LP batch solve workspace");
  }

  for(int system = 0; system < nsystem; system++) {
    NPLPSolvePolicyDiagnostics solve_diagnostics = {0, 0.0};
    size_t packed_index = 0;

    for(int a = 0; a < p; a++) {
      for(int b = a; b < p; b++) {
        const double value = packed_ptr[
          (size_t)system + (size_t)nsystem*packed_index
        ];
        workspace.gram_source[(size_t)a + (size_t)b*(size_t)p] = value;
        workspace.gram_source[(size_t)b + (size_t)a*(size_t)p] = value;
        packed_index++;
      }
    }

    if(TYPEOF(moments) == REALSXP) {
      const double *moment_matrix = REAL(moments);
      for(int term = 0; term < p; term++)
        workspace.rhs_source[term] = moment_matrix[
          (size_t)system + (size_t)nsystem*(size_t)term
        ];
    } else {
      for(int rhs = 0; rhs < nrhs; rhs++)
        for(int term = 0; term < p; term++)
        workspace.rhs_source[(size_t)term + (size_t)rhs*(size_t)p] = moment_ptr[term][
            (size_t)system + (size_t)nsystem*(size_t)rhs
          ];
    }

    status = np_lp_solve_workspace_solve_response(
      &workspace, p, nrhs, 1.0/mass_ptr[system], &solve_diagnostics);
    if(want_diagnostics) {
      INTEGER(ridge_steps)[system] = solve_diagnostics.ridge_steps;
      REAL(ridge_total)[system] = solve_diagnostics.ridge_total;
    }
    if(status != NP_LP_SOLVE_POLICY_OK) {
      failed_system = system + 1;
      break;
    }

    for(int rhs = 0; rhs < nrhs; rhs++) {
      double projected = 0.0;
      for(int term = 0; term < p; term++)
        projected += projection_ptr[term]*workspace.rhs_work[
          (size_t)term + (size_t)rhs*(size_t)p
        ];
      if(!R_FINITE(projected)) {
        status = NP_LP_SOLVE_POLICY_FINAL_FAILED;
        failed_system = system + 1;
        break;
      }
      value_ptr[(size_t)system + (size_t)rhs*(size_t)nsystem] = projected;
    }
    if(status != NP_LP_SOLVE_POLICY_OK)
      break;
  }

  np_lp_solve_workspace_clear(&workspace);

  result = PROTECT(np_lp_batch_result(
    (status == NP_LP_SOLVE_POLICY_OK) ? values : R_NilValue,
    status,
    failed_system,
    ridge_steps,
    ridge_total));
  nprotect++;
  UNPROTECT(nprotect);
  return result;
}
