.np_condens_engine_spec <- function(bws) {
  npConditionalRegEngineSpec(
    bws,
    where = "conditional density slice"
  )
}

.np_condens_is_already_proper_by_design <- function(bws) {
  .np_conditional_proper_certificate(bws, cdf = FALSE)
}

.np_conditional_proper_certificate <- function(bws, cdf) {
  spec <- .np_condens_engine_spec(bws)
  if (!(identical(spec$reg.engine, "lc") ||
        (identical(spec$reg.engine, "lp") && all(spec$degree.engine == 0L))))
    return(FALSE)
  nonnegative <- function(kernel, order) {
    identical(kernel, "uniform") ||
      (kernel %in% c("gaussian", "epanechnikov", "beta") && order == 2L)
  }
  # At fixed X, LC/LP0 must form a convex mixture. Higher-order predictor
  # kernels can give signed mixture weights even with a positive Y kernel.
  if (bws$xncon > 0L && !nonnegative(bws$cxkertype, bws$cxkerorder))
    return(FALSE)
  if (bws$yncon > 0L) {
    if (!nonnegative(bws$cykertype, bws$cykerorder) ||
        identical(bws$type, "generalized_nn"))
      return(FALSE)
    # Target-centred beta/boundary densities are not normalized as functions
    # of evaluation Y. Their observation-centred CDF siblings are proper.
    if (!cdf && (identical(bws$cykertype, "beta") ||
                 any(is.finite(bws$cykerlb[bws$iycon])) ||
                 any(is.finite(bws$cykerub[bws$iycon]))))
      return(FALSE)
  }
  # Both unordered density kernels are normalized over the retained support.
  if (bws$ynuno > 0L && !(bws$uykertype %in% c("aitchisonaitken", "liracine")) &&
      any(bws$ybw[bws$iyuno] != 0))
    return(FALSE)
  if (bws$ynord > 0L && !identical(bws$oykertype, "racineliyan") &&
      any(bws$ybw[bws$iyord] != 0))
    return(FALSE)
  TRUE
}

.np_condens_slice_dispatch_enabled <- function() {
  !identical(getOption("np.condens.proper.slice.enable"), FALSE)
}

.np_condens_validate_nonnegative_finite_numeric <- function(value, argname) {
  value <- as.double(value)[1L]
  if (!is.finite(value) || value < 0)
    stop(sprintf("'%s' must be a non-negative finite numeric scalar", argname))
  value
}
