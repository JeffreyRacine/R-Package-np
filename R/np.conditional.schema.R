# Plot grids must be constructed from the fitted conditional training schema.
# Preserve the conditional estimators' type policy, including categorical Y.
.np_conditional_training_types <- function(bws, xdat, ydat) {
  for (axis in c("x", "y")) {
    dat <- if (axis == "x") xdat else ydat
    argument <- paste0(axis, "dat")
    if (length(bws[[paste0(axis, "bw"), exact = TRUE]]) != ncol(dat))
      stop(sprintf("length of bandwidth vector does not match number of columns of '%s'",
                   argument), call. = FALSE)
    for (kind in c("con", "ord", "uno")) {
      mask <- bws[[paste0("i", axis, kind), exact = TRUE]]
      expected <- switch(kind, con = c("integer", "numeric"),
                          ord = "ordered", uno = "factor")
      if (any(mask) && !all(vapply(dat[, mask, drop = FALSE], inherits,
                                   logical(1), expected)))
        stop(sprintf("supplied bandwidths do not match '%s' in type", argument),
             call. = FALSE)
    }
  }
  invisible(NULL)
}
