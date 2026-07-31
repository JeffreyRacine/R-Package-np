#include <limits.h>
#include <stddef.h>

#include <R.h>
#include <Rinternals.h>

#include "categorical_profile_tile.h"

static int np_profile_matrix_dimensions(SEXP matrix,
                                        int *nrow,
                                        int *ncol)
{
  SEXP dim;

  if((nrow == NULL) || (ncol == NULL) ||
     (TYPEOF(matrix) != REALSXP) || !Rf_isMatrix(matrix))
    return 0;
  dim = Rf_getAttrib(matrix, R_DimSymbol);
  if((TYPEOF(dim) != INTSXP) || (XLENGTH(dim) != 2))
    return 0;
  *nrow = INTEGER(dim)[0];
  *ncol = INTEGER(dim)[1];
  return (*nrow > 0) && (*ncol >= 0);
}

static int np_profile_integer_vector_length(SEXP value, int expected)
{
  return (TYPEOF(value) == INTSXP) &&
    (expected >= 0) &&
    (XLENGTH(value) == (R_xlen_t)expected);
}

SEXP C_np_categorical_profile_kernel_tile(
  SEXP train_unordered,
  SEXP train_ordered,
  SEXP eval_unordered,
  SEXP eval_ordered,
  SEXP kernel_unordered,
  SEXP kernel_ordered,
  SEXP operator_code,
  SEXP lambda,
  SEXP num_categories,
  SEXP category_values,
  SEXP eval_start,
  SEXP eval_count)
{
  NPCategoricalProfileKernelSpec spec;
  NPCategoricalProfileTileStatus status;
  int ntrain_u, ntrain_o, neval_u, neval_o;
  int nunordered_train, nunordered_eval;
  int nordered_train, nordered_eval;
  int nvar;
  int start;
  int count;
  int i;
  size_t elements;
  size_t bytes;
  double **train_u_ptr;
  double **train_o_ptr;
  double **eval_u_ptr;
  double **eval_o_ptr;
  double **category_ptr;
  SEXP output;
  SEXP dim;

  if(!np_profile_matrix_dimensions(train_unordered,
                                   &ntrain_u,
                                   &nunordered_train) ||
     !np_profile_matrix_dimensions(train_ordered,
                                   &ntrain_o,
                                   &nordered_train) ||
     !np_profile_matrix_dimensions(eval_unordered,
                                   &neval_u,
                                   &nunordered_eval) ||
     !np_profile_matrix_dimensions(eval_ordered,
                                   &neval_o,
                                   &nordered_eval))
    Rf_error("categorical-profile inputs must be numeric matrices");
  if((ntrain_u != ntrain_o) || (neval_u != neval_o) ||
     (nunordered_train != nunordered_eval) ||
     (nordered_train != nordered_eval))
    Rf_error("categorical-profile train/evaluation dimensions are inconsistent");
  if(nunordered_train > INT_MAX - nordered_train)
    Rf_error("categorical-profile variable count exceeds native capacity");
  nvar = nunordered_train + nordered_train;
  if(nvar <= 0)
    Rf_error("categorical-profile inputs require at least one variable");

  if(!np_profile_integer_vector_length(kernel_unordered,
                                       nunordered_train) ||
     !np_profile_integer_vector_length(kernel_ordered,
                                       nordered_train) ||
     !np_profile_integer_vector_length(operator_code, nvar) ||
     !np_profile_integer_vector_length(num_categories, nvar))
    Rf_error("categorical-profile kernel metadata has inconsistent length");
  if((TYPEOF(lambda) != REALSXP) ||
     (XLENGTH(lambda) != (R_xlen_t)nvar))
    Rf_error("categorical-profile bandwidths must be a numeric vector of variable length");
  if((TYPEOF(category_values) != VECSXP) ||
     (XLENGTH(category_values) != (R_xlen_t)nvar))
    Rf_error("categorical-profile supports must be a list of variable length");
  if((XLENGTH(eval_start) != 1) || (XLENGTH(eval_count) != 1))
    Rf_error("categorical-profile tile range must use scalar integers");
  start = Rf_asInteger(eval_start);
  count = Rf_asInteger(eval_count);
  if((start == NA_INTEGER) || (count == NA_INTEGER))
    Rf_error("categorical-profile tile range must use finite scalar integers");
  if(start <= 0)
    Rf_error("categorical-profile evaluation start must be positive");
  start--;
  if((count <= 0) || (start > neval_u) ||
     (count > neval_u - start))
    Rf_error("categorical-profile evaluation range is invalid");

  train_u_ptr = (double **)R_alloc(
    (size_t)(nunordered_train > 0 ? nunordered_train : 1),
                                   sizeof(*train_u_ptr));
  train_o_ptr = (double **)R_alloc(
    (size_t)(nordered_train > 0 ? nordered_train : 1),
                                   sizeof(*train_o_ptr));
  eval_u_ptr = (double **)R_alloc(
    (size_t)(nunordered_train > 0 ? nunordered_train : 1),
                                  sizeof(*eval_u_ptr));
  eval_o_ptr = (double **)R_alloc(
    (size_t)(nordered_train > 0 ? nordered_train : 1),
                                  sizeof(*eval_o_ptr));
  category_ptr = (double **)R_alloc((size_t)nvar,
                                    sizeof(*category_ptr));

  for(i = 0; i < nunordered_train; i++) {
    train_u_ptr[i] = REAL(train_unordered) + (size_t)i*(size_t)ntrain_u;
    eval_u_ptr[i] = REAL(eval_unordered) + (size_t)i*(size_t)neval_u;
  }
  for(i = 0; i < nordered_train; i++) {
    train_o_ptr[i] = REAL(train_ordered) + (size_t)i*(size_t)ntrain_o;
    eval_o_ptr[i] = REAL(eval_ordered) + (size_t)i*(size_t)neval_o;
  }
  for(i = 0; i < nvar; i++) {
    SEXP support = VECTOR_ELT(category_values, i);
    const int ncat = INTEGER(num_categories)[i];

    category_ptr[i] = NULL;
    if(ncat <= 0)
      Rf_error("categorical-profile category counts must be positive");
    if(i < nunordered_train) {
      if(support != R_NilValue) {
        if((TYPEOF(support) != REALSXP) ||
           (XLENGTH(support) != (R_xlen_t)ncat))
          Rf_error("unordered categorical-profile support has inconsistent length");
        category_ptr[i] = REAL(support);
      }
    } else {
      if((TYPEOF(support) != REALSXP) ||
         (XLENGTH(support) != (R_xlen_t)ncat))
        Rf_error("ordered categorical-profile support has inconsistent length");
      category_ptr[i] = REAL(support);
    }
  }

  status = np_categorical_profile_tile_bytes((size_t)ntrain_u,
                                             (size_t)count,
                                             &elements,
                                             &bytes);
  if(status != NP_PROFILE_TILE_OK)
    Rf_error("%s", np_categorical_profile_tile_status_message(status));
  if(bytes > NP_CATEGORICAL_PROFILE_TILE_MAX_BYTES)
    Rf_error("categorical-profile tile exceeds the 64 MiB output ceiling");
  if(elements > (size_t)R_XLEN_T_MAX)
    Rf_error("categorical-profile tile exceeds the R vector-length limit");

  output = PROTECT(Rf_allocVector(REALSXP, (R_xlen_t)elements));
  spec.ntrain = ntrain_u;
  spec.neval = neval_u;
  spec.nunordered = nunordered_train;
  spec.nordered = nordered_train;
  spec.train_unordered = train_u_ptr;
  spec.train_ordered = train_o_ptr;
  spec.eval_unordered = eval_u_ptr;
  spec.eval_ordered = eval_o_ptr;
  spec.kernel_unordered = INTEGER(kernel_unordered);
  spec.kernel_ordered = INTEGER(kernel_ordered);
  spec.operator_code = INTEGER(operator_code);
  spec.lambda = REAL(lambda);
  spec.num_categories = INTEGER(num_categories);
  spec.category_values = category_ptr;

  status = np_categorical_profile_tile_fill(&spec,
                                            start,
                                            count,
                                            REAL(output),
                                            elements);
  if(status != NP_PROFILE_TILE_OK) {
    UNPROTECT(1);
    Rf_error("%s", np_categorical_profile_tile_status_message(status));
  }

  dim = PROTECT(Rf_allocVector(INTSXP, 2));
  INTEGER(dim)[0] = ntrain_u;
  INTEGER(dim)[1] = count;
  Rf_setAttrib(output, R_DimSymbol, dim);
  UNPROTECT(2);
  return output;
}
