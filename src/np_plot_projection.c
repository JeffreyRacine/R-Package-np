/* Plot-only contraction of an accepted fit projection and kernel row.
 * No fit, rank decision or factorization occurs in this helper. */
#include <R.h>
#include <Rinternals.h>
SEXP C_np_plot_projection_apply(SEXP basis, SEXP projection, SEXP kernel)
{
  if(TYPEOF(basis) != REALSXP || !isMatrix(basis) ||
     TYPEOF(projection) != REALSXP || !isMatrix(projection))
    error("plot projection requires numeric matrices");
  const int n = nrows(basis), p = ncols(basis), g = nrows(projection);
  if(n < 1 || p < 1 || g < 1 || ncols(projection) != p ||
     (kernel != R_NilValue && (TYPEOF(kernel) != REALSXP || !isMatrix(kernel) ||
                              nrows(kernel) != n || ncols(kernel) != g)))
    error("plot projection dimensions do not agree");
  SEXP out = PROTECT(allocMatrix(REALSXP, g, n));
  for(int j = 0; j < g; ++j) {
    R_CheckUserInterrupt();
    for(int i = 0; i < n; ++i) {
      double value = 0.0;
      for(int k = 0; k < p; ++k)
        value += REAL(basis)[i + (size_t)k*n] *
                 REAL(projection)[j + (size_t)k*g];
      if(kernel != R_NilValue) value *= REAL(kernel)[i + (size_t)j*n];
      if(!R_FINITE(value)) error("plot projection has nonfinite influence");
      REAL(out)[j + (size_t)i*g] = value;
    }
  }
  UNPROTECT(1);
  return out;
}
