#line 11113 "np.c"
SEXP C_np_regression_lp_apply_conditional_ctx(SEXP txuno,
                                              SEXP txord,
                                              SEXP txcon,
                                              SEXP exuno,
                                              SEXP exord,
                                              SEXP excon,
                                              SEXP rhs,
                                              SEXP rbw,
                                              SEXP bwtype,
                                              SEXP kernel_x,
                                              SEXP kernel_xu,
                                              SEXP kernel_xo,
                                              SEXP use_tree,
                                              SEXP glp_degree,
                                              SEXP glp_gradient_order,
                                              SEXP glp_bernstein,
                                              SEXP glp_basis,
                                              SEXP continuous_kernel_family,
                                              SEXP continuous_kernel_order,
                                              SEXP ckerlb,
                                              SEXP ckerub,
                                              SEXP categorical_compress,
                                              SEXP return_hat,
                                              SEXP train_is_eval,
                                              SEXP leave_one_out,
                                              SEXP allow_empty_rows,
                                              SEXP declared_support)
{
#ifdef NP_CF316_LITERAL_ENTRIES
  if(asLogical(return_hat) != TRUE || asLogical(leave_one_out) != FALSE)
    error("conditional literal replay requires a full hat and no leave-one-out request");
#endif

  return np_regression_lp_apply_conditional_impl(
    txuno, txord, txcon, exuno, exord, excon, rhs, rbw, bwtype,
    kernel_x, kernel_xu, kernel_xo, use_tree, glp_degree,
    glp_gradient_order, glp_bernstein, glp_basis,
    continuous_kernel_family, continuous_kernel_order, ckerlb, ckerub,
    categorical_compress, return_hat, train_is_eval, leave_one_out,
    R_NilValue, R_NilValue, R_NilValue, R_NilValue, R_NilValue,
    R_NilValue, allow_empty_rows, declared_support);
}
