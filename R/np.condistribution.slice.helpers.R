.np_condist_engine_spec <- function(bws) {
  npConditionalRegEngineSpec(
    bws,
    where = "conditional distribution slice"
  )
}

.np_condist_is_already_proper_by_design <- function(bws) {
  spec <- .np_condist_engine_spec(bws)

  identical(spec$reg.engine, "lc") ||
    (identical(spec$reg.engine, "lp") && all(spec$degree.engine == 0L))
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
