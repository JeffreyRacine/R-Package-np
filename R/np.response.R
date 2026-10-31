# Validate the effective public outcome without coercing, copying or scanning it.
# Predictor scores and operator right-hand sides have separate internal owners.
.np_require_numeric_response <- function(y, where, argument) {
  if (!is.numeric(y) || !(is.integer(y) || is.double(y)) ||
      is.factor(y) || !is.null(dim(y)))
    stop(paste0(where, "(): '", argument,
                "' must be a numeric (integer or double) response vector; ",
                "factors and ordered factors are not supported."), call. = FALSE)
  invisible(NULL)
}
