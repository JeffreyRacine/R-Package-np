.np_condist_engine_spec <- function(bws) {
  reg.engine <- if (is.null(bws$regtype.engine)) {
    if (is.null(bws$regtype)) "lc" else as.character(bws$regtype)
  } else {
    as.character(bws$regtype.engine)
  }

  degree.engine <- if (is.null(bws$degree.engine)) {
    if (isTRUE(bws$xncon > 0L)) {
      if (identical(reg.engine, "lc")) {
        rep.int(0L, bws$xncon)
      } else {
        as.integer(npValidateGlpDegree(
          regtype = "lp",
          degree = bws$degree,
          ncon = bws$xncon
        ))
      }
    } else {
      integer(0L)
    }
  } else {
    as.integer(bws$degree.engine)
  }

  list(
    reg.engine = reg.engine,
    degree.engine = degree.engine
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
