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
    if(NP_CF316_DISC_UPPER(code, c->lambda[j], c->counts[j])) {
      const double a = ukf[code](1, c->lambda[j], c->counts[j]);
      const double b = ukf[code](0, c->lambda[j], c->counts[j]);
      if(NP_CF316_DISC_CONSTANT(a, b)) {
        c->uno_ok[j] = 1; c->uno_const[j] = 0.5*(a+b);
      }
    }
  }
  for(int j = 0; j < no; ++j) {
    const int code = c->codes[nc+nu+j], index = nu+j;
    if(code < 0 || code >= 16) error("plot gate ordered kernel code is invalid");
    if(NP_CF316_DISC_ORDERED_UPPER(code, c->lambda[index])) {
      double a, b;
      np_ordered_kernel_range(code, c->lambda[index], c->support[index],
                              c->counts[index], &a, &b);
      if(NP_CF316_DISC_CONSTANT(a, b)) {
        c->ord_ok[j] = 1; c->ord_const[j] = 0.5*(a+b);
      }
    }
  }
  np_gate_ctx_set(&c->gate, nc, nu, no, c->codes, c->codes+nc,
    c->codes+nc+nu, c->op, c->cont_ok, c->cont_hmin, c->cont_k0,
    c->uno_ok, c->uno_const, c->ord_ok, c->ord_const);
  if(!np_gate_ctx_is_sane(&c->gate)) error("plot gate prepared context is invalid");
}
