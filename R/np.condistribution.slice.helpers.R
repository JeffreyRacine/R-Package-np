.np_condist_engine_spec <- function(bws) {
  npConditionalRegEngineSpec(
    bws,
    where = "conditional distribution slice"
  )
}

.np_condist_is_already_proper_by_design <- function(bws) {
  .np_conditional_proper_certificate(bws, cdf = TRUE)
}

.np_condist_slice_dispatch_enabled <- function() {
  !identical(getOption("np.condist.proper.slice.enable"), FALSE)
}

.np_condist_validate_nonnegative_finite_numeric <- function(value, argname) {
  value <- as.double(value)[1L]
  if (!is.finite(value) || value < 0)
    stop(sprintf("'%s' must be a non-negative finite numeric scalar", argname))
  value
}
