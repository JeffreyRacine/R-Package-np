# The native interface consumes explicit vectors/matrices/data frames. Inspect
# only argument names: an unsupported data promise (including NULL) is not input.
.np_native_data_error <- function() {
  stop("`data=` is supported only with the formula interface. For native calls, supply vectors, matrices or data frames directly, or use `with(data, ...)`.",
       call. = FALSE)
}

.np_reject_native_data <- function(dots) {
  if ("data" %in% names(dots)) .np_native_data_error()
  invisible(NULL)
}

# A bare unbound dispatch symbol cannot select any valid S3 method. Improve the
# common IV misuse without forcing y/data, inspecting a data mask, or catching
# arbitrary dispatch errors. Bound formula symbols and factories use UseMethod.
.np_iv_native_data_preflight <- function(y.expr, dots, caller) {
  if ("data" %in% names(dots) && all(c("z", "w") %in% names(dots)) &&
      is.symbol(y.expr) && nzchar(as.character(y.expr)) &&
      !exists(as.character(y.expr), envir = caller, inherits = TRUE))
    .np_native_data_error()
  invisible(NULL)
}

.np_plot_data_error <- function() {
  stop("This plot method does not consume `data=` on the selected route. Refit the model with the intended data and plot that model, or supply the supported explicit plotting inputs.",
       call. = FALSE)
}

# match.call drops explicitly empty actuals when their formals have defaults.
# Replace only those syntax slots before matching; never evaluate an argument.
# Ellipsis in the inspected call belongs to its caller, not this helper caller.
.np_data_formal_supplied <- function(call, definition) {
  actuals <- as.list(call)
  empty <- vapply(actuals, identical, logical(1L), quote(expr = ))
  actuals[empty] <- list(quote(.np_missing_actual))
  "data" %in% names(match.call(definition = definition,
                             call = as.call(actuals), expand.dots = FALSE,
                             envir = parent.frame(2L)))
}
