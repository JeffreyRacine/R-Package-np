/* Private plot-inference dispatcher. Include in jksum.c after its gate/kernel helpers.
 * No existing function or protected regression owner is modified.
 * Physical bandwidth and complete declared support are supplied by the
 * existing R decoders. This wrapper owns only its new-request gate state.
 */
extern SEXP C_np_kernelsum(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
  SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
  SEXP, SEXP, SEXP);
extern SEXP C_np_regression_lp_apply_conditional_ctx(SEXP, SEXP, SEXP,
  SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
  SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP,
  SEXP, SEXP, SEXP, SEXP);

typedef struct {
  SEXP args, geometry;
  int dispatch, installed, n, m, dims[3], child_n, child_m;
  int base[3], bandwidth_type, support_columns;
  const int *parsed_counts;
  const double *const *parsed_support;
  NP_GateOverrideCtx saved, gate;
  int *codes, *op, *cont_ok, *uno_ok, *ord_ok, *counts;
  double *cont_hmin, *cont_k0, *uno_const, *ord_const, *lambda;
  double **train_con, **eval_con, **support;
} NPPlotGateCall;

static void *np_plot_gate_alloc(size_t count, size_t size)
{
  if(count == 0) count = 1;
  if(size == 0 || count > SIZE_MAX / size)
    error("plot gate scope allocation overflow");
  void *p = calloc(count, size);
  if(p == NULL) error("plot gate scope allocation failed");
  return p;
}

/* Called for success and non-local exit; restore before freeing owned arrays.
 * The displaced context may borrow an earlier scope's arrays: do not free it.
 */
static void np_plot_gate_cleanup(void *data, Rboolean jump)
{
  NPPlotGateCall *c = (NPPlotGateCall *)data;
  (void)jump;
  if(c->installed) {
    np_gate_override_ctx = c->saved;
    c->installed = 0;
  }
  if(c->support != NULL)
    for(int j = 0; j < c->dims[0] + c->dims[1]; ++j) free(c->support[j]);
  free(c->support); free(c->counts); free(c->lambda);
  free(c->codes); free(c->op); free(c->cont_ok);
  free(c->cont_hmin); free(c->cont_k0);
  free(c->uno_ok); free(c->uno_const); free(c->ord_ok); free(c->ord_const);
  free(c->train_con); free(c->eval_con);
}

static void np_plot_gate_names(SEXP x, const char *const *names, int n,
                               int required)
{
  if(TYPEOF(x) != VECSXP || XLENGTH(x) != n)
    error("plot gate scope list shape mismatch");
  SEXP actual = getAttrib(x, R_NamesSymbol);
  if(!required && actual == R_NilValue) return;
  if(TYPEOF(actual) != STRSXP || XLENGTH(actual) != n)
    error("plot gate scope requires exact field names");
  for(int j = 0; j < n; ++j)
    if(STRING_ELT(actual, j) == NA_STRING ||
       strcmp(CHAR(STRING_ELT(actual, j)), names[j]) != 0)
      error("plot gate scope field order mismatch");
}

static int np_plot_gate_integer(SEXP x)
{
  if(TYPEOF(x) != INTSXP || XLENGTH(x) != 1 || INTEGER(x)[0] == NA_INTEGER)
    error("plot gate scope requires an integer scalar");
  return INTEGER(x)[0];
}

static void np_plot_gate_matrix(SEXP x, int *n, int *p)
{
  SEXP dim = getAttrib(x, R_DimSymbol);
  if(TYPEOF(x) != REALSXP || TYPEOF(dim) != INTSXP || XLENGTH(dim) != 2 ||
     INTEGER(dim)[0] <= 0 || INTEGER(dim)[1] < 0)
    error("plot gate geometry requires double matrices, including zero-column groups");
  *n = INTEGER(dim)[0]; *p = INTEGER(dim)[1];
  if(XLENGTH(x) != (R_xlen_t)*n * *p)
    error("plot gate geometry matrix extent mismatch");
  for(R_xlen_t i = 0; i < XLENGTH(x); ++i)
    if(!R_FINITE(REAL(x)[i])) error("plot gate geometry must be finite");
}

static void np_plot_gate_index(SEXP index, int n, int limit)
{
  if(TYPEOF(index) != INTSXP || XLENGTH(index) != n)
    error("plot gate geometry index shape mismatch");
  for(int i = 0; i < n; ++i)
    if(INTEGER(index)[i] < 1 || INTEGER(index)[i] > limit ||
       INTEGER(index)[i] == NA_INTEGER)
      error("plot gate geometry index is out of range");
}

static void np_plot_gate_mapped_values(SEXP child, int n, int p,
                                       SEXP complete, int full_n, SEXP index)
{
  if(TYPEOF(child) != REALSXP ||
     (p > 0 && XLENGTH(child) != (R_xlen_t)n*p) ||
     (p == 0 && XLENGTH(child) > 1))
    error("plot gate child numeric buffer shape mismatch");
  for(int j = 0; j < p; ++j)
    for(int i = 0; i < n; ++i)
      if(!R_FINITE(REAL(child)[(R_xlen_t)j*n+i]) ||
         REAL(child)[(R_xlen_t)j*n+i] !=
         REAL(complete)[(R_xlen_t)j*full_n+INTEGER(index)[i]-1])
        error("plot gate child rows do not match complete geometry");
}

/* np.c supplies support parsed by its existing canonical validator. */
static void np_plot_gate_support(NPPlotGateCall *c)
{
  const int nd = c->dims[0]+c->dims[1];
  if(c->support_columns != nd ||
     (nd > 0 && (c->parsed_counts == NULL || c->parsed_support == NULL)))
    error("plot gate parsed support shape mismatch");
  c->counts = np_plot_gate_alloc((size_t)nd, sizeof(int));
  c->support = np_plot_gate_alloc((size_t)nd, sizeof(double *));
  for(int j = 0; j < nd; ++j) {
    if(c->parsed_counts[j] < 1 || c->parsed_support[j] == NULL)
      error("plot gate parsed support is incomplete");
    c->counts[j] = c->parsed_counts[j];
    c->support[j] = np_plot_gate_alloc((size_t)c->counts[j], sizeof(double));
    memcpy(c->support[j], c->parsed_support[j], (size_t)c->counts[j]*sizeof(double));
  }
}

#include "conditional_literal/jksum_np_plot_gate_prepare.inc"

static void np_plot_gate_prepare_conditional_literal(NPPlotGateCall *c);

#include "conditional_literal/jksum_np_plot_gate_body.inc"

static SEXP np_plot_gate_body_conditional_literal(void *data);

#include "conditional_literal/jksum_np_plot_gate_dispatch_native.inc"

SEXP np_plot_gate_dispatch_native_conditional_literal(SEXP kind, SEXP native_args, SEXP geometry,
  int support_columns, const int *support_counts,
  const double *const *support_values);
