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

static void np_plot_gate_prepare(NPPlotGateCall *c)
{
  static const char *const geometry_names[] = {"train", "eval",
    "physical.bandwidth", "support", "train.index", "eval.index"};
  static const char *const group_names[] = {"uno", "ord", "con"};
  static const char *const sum_names[] = {"tuno", "tord", "tcon", "ty",
    "weights", "euno", "eord", "econ", "bw", "mcv", "padnum", "op",
    "myopti", "kpow", "ksum_len", "pksum_len", "kw_len", "ckerlb", "ckerub"};
  static const char *const lp_names[] = {"txuno", "txord", "txcon", "exuno",
    "exord", "excon", "rhs", "rbw", "bwtype", "kernel_x", "kernel_xu",
    "kernel_xo", "use_tree", "glp_degree", "glp_gradient_order", "glp_bernstein",
    "glp_basis", "continuous_kernel_family", "continuous_kernel_order", "ckerlb",
    "ckerub", "categorical_compress", "return_hat", "train_is_eval",
    "leave_one_out", "allow_empty_rows", "declared_support"};
  np_plot_gate_names(c->geometry, geometry_names, 6, 1);
  np_plot_gate_names(c->args, c->dispatch ? lp_names : sum_names,
                      c->dispatch ? 27 : 19, 0);
  SEXP train = VECTOR_ELT(c->geometry, 0), eval = VECTOR_ELT(c->geometry, 1);
  SEXP bw = VECTOR_ELT(c->geometry, 2), support = VECTOR_ELT(c->geometry, 3);
  SEXP ti = VECTOR_ELT(c->geometry, 4), ei = VECTOR_ELT(c->geometry, 5);
  np_plot_gate_names(train, group_names, 3, 1);
  np_plot_gate_names(eval, group_names, 3, 1);
  for(int group = 0; group < 3; ++group) {
    int n, m, p, q;
    np_plot_gate_matrix(VECTOR_ELT(train, group), &n, &p);
    np_plot_gate_matrix(VECTOR_ELT(eval, group), &m, &q);
    if(p != q || (group > 0 && (n != c->n || m != c->m)))
      error("plot gate complete geometry dimensions differ");
    c->n = n; c->m = m; c->dims[group] = p;
  }
  const int nu = c->dims[0], no = c->dims[1], nc = c->dims[2];
  const R_xlen_t total = (R_xlen_t)nu+no+nc;
  if(total <= 0 || total > INT_MAX || TYPEOF(bw) != REALSXP || XLENGTH(bw) != total)
    error("plot gate physical bandwidth dimensions differ");
  c->codes = np_plot_gate_alloc((size_t)total, sizeof(int));
  c->op = np_plot_gate_alloc((size_t)total, sizeof(int));
  c->lambda = np_plot_gate_alloc((size_t)(nu+no), sizeof(double));
  if(TYPEOF(support) != VECSXP || XLENGTH(support) != nu+no)
    error("plot gate declared support shape mismatch");
  for(int j = 0; j < XLENGTH(c->args); ++j) {
    const int type = TYPEOF(VECTOR_ELT(c->args, j));
    if(c->dispatch && j == 26) {
      if(type != VECSXP) error("plot gate LP requires a decoded support list");
    } else if(type != REALSXP && type != INTSXP && type != LGLSXP)
      error("plot gate native payload must contain ordinary numeric vectors");
  }
  SEXP child_bw;
  int identity = 0;
  if(c->dispatch == 0) {
    SEXP options = VECTOR_ELT(c->args, 12), op = VECTOR_ELT(c->args, 11);
    if(TYPEOF(options) != INTSXP || XLENGTH(options) < KWS_OPTIONS_COUNT ||
       TYPEOF(op) != INTSXP || XLENGTH(op) != total)
      error("plot gate kernel-sum option shape mismatch");
    const int *o = INTEGER(options);
    if(o[KWS_NUNOI] != nu || o[KWS_NORDI] != no || o[KWS_NCONI] != nc ||
       o[KWS_LSFI] != SF_ARB || o[KWS_TNOBSI] <= 0 || o[KWS_ENOBSI] <= 0 ||
       (o[KWS_TISEI] != 0 && o[KWS_TISEI] != 1))
      error("plot gate kernel-sum requires matching physical bandwidth marshalling");
    c->child_n = o[KWS_TNOBSI]; c->child_m = o[KWS_ENOBSI];
    c->base[0] = o[KWS_UKRNEVI]; c->base[1] = o[KWS_OKRNEVI];
    c->base[2] = o[KWS_CKRNEVI]; c->bandwidth_type = o[KWS_BWI];
    identity = o[KWS_TISEI]; child_bw = VECTOR_ELT(c->args, 8);
    for(int j = 0; j < total; ++j) {
      if(INTEGER(op)[j] != OP_NORMAL)
        error("plot gate scratch kernel-sum entry currently requires normal operators");
      c->op[j] = INTEGER(op)[j];
    }
  } else {
    for(int group = 0; group < 3; ++group) {
      int n, m, p, q;
      np_plot_gate_matrix(VECTOR_ELT(c->args, group), &n, &p);
      np_plot_gate_matrix(VECTOR_ELT(c->args, group+3), &m, &q);
      if(p != c->dims[group] || q != p ||
         (group > 0 && (n != c->child_n || m != c->child_m)))
        error("plot gate LP child matrix dimensions differ");
      c->child_n = n; c->child_m = m;
    }
    c->base[2] = np_plot_gate_integer(VECTOR_ELT(c->args, 9));
    c->base[0] = np_plot_gate_integer(VECTOR_ELT(c->args, 10));
    c->base[1] = np_plot_gate_integer(VECTOR_ELT(c->args, 11));
    c->bandwidth_type = np_plot_gate_integer(VECTOR_ELT(c->args, 8));
    child_bw = VECTOR_ELT(c->args, 7);
    /* calloc initialized the LP owner's main operator to OP_NORMAL. */
    SEXP child_support = VECTOR_ELT(c->args, 26);
    if(TYPEOF(child_support) != VECSXP || XLENGTH(child_support) != nu+no)
      error("plot gate LP support differs from complete geometry");
    for(int j = 0; j < nu+no; ++j) {
      SEXP a = VECTOR_ELT(child_support, j), b = VECTOR_ELT(support, j);
      if(TYPEOF(a) != REALSXP || TYPEOF(b) != REALSXP || XLENGTH(a) != XLENGTH(b) ||
         memcmp(REAL(a), REAL(b), (size_t)XLENGTH(a)*sizeof(double)) != 0)
        error("plot gate LP decoded support differs from complete geometry");
    }
  }
  if(c->bandwidth_type != BW_FIXED && c->bandwidth_type != BW_GEN_NN &&
     c->bandwidth_type != BW_ADAP_NN)
    error("plot gate bandwidth type is invalid");
  if(TYPEOF(child_bw) != REALSXP || XLENGTH(child_bw) != total)
    error("plot gate child bandwidth shape mismatch");
  for(int j = 0; j < total; ++j) {
    if(!R_FINITE(REAL(bw)[j]) || REAL(bw)[j] != REAL(child_bw)[j])
      error("plot gate physical bandwidth does not match native payload");
    if(j >= nc) c->lambda[j-nc] = REAL(bw)[j];
  }
  np_plot_gate_index(ti, c->child_n, c->n);
  np_plot_gate_index(ei, c->child_m, c->m);
  /* This request tiles evaluation only. A donor subset/permutation changes
   * the private NN geometry contract; do not admit it implicitly. A caller
   * proving donor-order invariance supplies a correspondingly permuted FULL
   * parent training frame and its original occurrence maps instead. */
  if(c->child_n != c->n)
    error("plot gate scope requires the complete training donor set");
  for(int i = 0; i < c->n; ++i)
    if(INTEGER(ti)[i] != i+1)
      error("plot gate scope requires identity training order");
  if(identity && c->child_n != c->child_m)
    error("plot gate implicit training identity has unequal row counts");
  for(int group = 0; group < 3; ++group) {
    SEXP child_train = VECTOR_ELT(c->args, group);
    SEXP child_eval = identity ? child_train :
      VECTOR_ELT(c->args, group+(c->dispatch ? 3 : 5));
    np_plot_gate_mapped_values(child_train, c->child_n, c->dims[group],
                               VECTOR_ELT(train, group), c->n, ti);
    np_plot_gate_mapped_values(child_eval, c->child_m, c->dims[group],
                               VECTOR_ELT(eval, group), c->m, ei);
  }
  np_plot_gate_support(c);
  if(c->dispatch == 0 && nu+no > 0) {
    const int rows = INTEGER(VECTOR_ELT(c->args, 12))[KWS_MLEVI];
    SEXP mcv = VECTOR_ELT(c->args, 9), pad = VECTOR_ELT(c->args, 10);
    if(rows <= 0 || TYPEOF(mcv) != REALSXP ||
       XLENGTH(mcv) != (R_xlen_t)rows*(nu+no) ||
       TYPEOF(pad) != REALSXP || XLENGTH(pad) != 1 || !R_FINITE(REAL(pad)[0]))
      error("plot gate kernel-sum declared support buffer differs");
    for(int j = 0; j < nu+no; ++j) {
      if(c->counts[j] > rows) error("plot gate support exceeds kernel-sum buffer");
      for(int k = 0; k < rows; ++k) {
        const double expected = k < c->counts[j] ? c->support[j][k] : REAL(pad)[0];
        if(REAL(mcv)[(R_xlen_t)j*rows+k] != expected)
          error("plot gate kernel-sum support is not the complete decoded support");
      }
    }
  }
  for(int j = 0; j < nc; ++j) c->codes[j] = c->base[2];
  for(int j = 0; j < nu; ++j) c->codes[nc+j] = c->base[0];
  for(int j = 0; j < no; ++j) c->codes[nc+nu+j] = c->base[1];

  /* Refresh the same option snapshot before installing any override. */
  np_refresh_runtime_tolerances();
  c->train_con = np_plot_gate_alloc((size_t)nc, sizeof(double *));
  c->eval_con = np_plot_gate_alloc((size_t)nc, sizeof(double *));
  for(int j = 0; j < nc; ++j) {
    c->train_con[j] = REAL(VECTOR_ELT(train, 2))+(R_xlen_t)j*c->n;
    c->eval_con[j] = REAL(VECTOR_ELT(eval, 2))+(R_xlen_t)j*c->m;
  }
  if(nc > 0) {
    if(np_largeh_enabled()) {
      if(!np_cont_largeh_build_params(c->n, c->m, nc, c->codes,
           c->train_con, c->eval_con, np_cont_largeh_rel_tol(),
           &c->cont_ok, &c->cont_hmin, &c->cont_k0))
        error("plot gate continuous geometry preparation failed");
    } else {
      c->cont_ok = np_plot_gate_alloc((size_t)nc, sizeof(int));
      c->cont_hmin = np_plot_gate_alloc((size_t)nc, sizeof(double));
      c->cont_k0 = np_plot_gate_alloc((size_t)nc, sizeof(double));
      for(int j = 0; j < nc; ++j) c->cont_hmin[j] = DBL_MAX;
    }
  }
  c->uno_ok = np_plot_gate_alloc((size_t)nu, sizeof(int));
  c->uno_const = np_plot_gate_alloc((size_t)nu, sizeof(double));
  c->ord_ok = np_plot_gate_alloc((size_t)no, sizeof(int));
  c->ord_const = np_plot_gate_alloc((size_t)no, sizeof(double));
  double (*const ukf[])(int, double, int) = {np_uaa, np_unli_racine,
    np_econvol_uaa, np_econvol_unli_racine, np_score_uaa, np_score_unli_racine};
  const int nuk = (int)(sizeof(ukf)/sizeof(ukf[0]));
  for(int j = 0; j < nu; ++j) {
    const int code = c->codes[nc+j];
    if(code < 0 || code >= nuk) error("plot gate unordered kernel code is invalid");
    if(np_disc_near_upper(code, c->lambda[j], c->counts[j])) {
      const double a = ukf[code](1, c->lambda[j], c->counts[j]);
      const double b = ukf[code](0, c->lambda[j], c->counts[j]);
      if(np_disc_near_const_kernel(a, b)) {
        c->uno_ok[j] = 1; c->uno_const[j] = 0.5*(a+b);
      }
    }
  }
  for(int j = 0; j < no; ++j) {
    const int code = c->codes[nc+nu+j], index = nu+j;
    if(code < 0 || code >= 16) error("plot gate ordered kernel code is invalid");
    if(np_disc_ordered_near_upper(code, c->lambda[index])) {
      double a, b;
      np_ordered_kernel_range(code, c->lambda[index], c->support[index],
                              c->counts[index], &a, &b);
      if(np_disc_near_const_kernel(a, b)) {
        c->ord_ok[j] = 1; c->ord_const[j] = 0.5*(a+b);
      }
    }
  }
  np_gate_ctx_set(&c->gate, nc, nu, no, c->codes, c->codes+nc,
    c->codes+nc+nu, c->op, c->cont_ok, c->cont_hmin, c->cont_k0,
    c->uno_ok, c->uno_const, c->ord_ok, c->ord_const);
  if(!np_gate_ctx_is_sane(&c->gate)) error("plot gate prepared context is invalid");
}

static SEXP np_plot_gate_body(void *data)
{
  NPPlotGateCall *c = (NPPlotGateCall *)data;
  np_plot_gate_prepare(c);
  c->saved = np_gate_override_ctx;
  c->installed = 1;
  np_gate_override_ctx = c->gate;
  /* Fixed native dispatch only: no eval(), R closure, or callback argument. */
#define A(i) VECTOR_ELT(c->args, (i))
  if(c->dispatch == 0)
    return C_np_kernelsum(A(0), A(1), A(2), A(3), A(4), A(5), A(6), A(7),
      A(8), A(9), A(10), A(11), A(12), A(13), A(14), A(15), A(16), A(17), A(18));
  return C_np_regression_lp_apply_conditional_ctx(A(0), A(1), A(2), A(3),
    A(4), A(5), A(6), A(7), A(8), A(9), A(10), A(11), A(12), A(13),
    A(14), A(15), A(16), A(17), A(18), A(19), A(20), A(21), A(22),
    A(23), A(24), A(25), A(26));
#undef A
}

SEXP np_plot_gate_dispatch_native(SEXP kind, SEXP native_args, SEXP geometry,
  int support_columns, const int *support_counts,
  const double *const *support_values)
{
  NPPlotGateCall call = {0};
  if(TYPEOF(kind) != STRSXP || XLENGTH(kind) != 1 || STRING_ELT(kind, 0) == NA_STRING)
    error("plot gate dispatcher requires one fixed native operation");
  const char *label = CHAR(STRING_ELT(kind, 0));
  if(strcmp(label, "kernelsum") == 0) call.dispatch = 0;
  else if(strcmp(label, "lp_apply_ctx") == 0) call.dispatch = 1;
  else error("unknown plot gate native operation");
  call.args = native_args; call.geometry = geometry;
  call.support_columns = support_columns; call.parsed_counts = support_counts;
  call.parsed_support = support_values;
  return R_UnwindProtect(np_plot_gate_body, &call,
                         np_plot_gate_cleanup, &call, NULL);
}
