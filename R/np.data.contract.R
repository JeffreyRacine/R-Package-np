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

# Call immediately before the incumbent owner first forces this dot. Only a
# directly supplied named symbol has an unambiguous lexical lookup here.
# Forwarded promises, factories and existing bindings keep ordinary evaluation.
.np_data_argument_preflight <- function(expr, argument, call, caller, value) {
  if (missing(expr) || !is.symbol(expr) || !nzchar(as.character(expr)) ||
      !nzchar(argument) || sum(names(call) == argument) != 1L ||
      !identical(call[[argument]], expr) ||
      exists(as.character(expr), envir = caller, inherits = TRUE))
    return(invisible(NULL))
  # Absence of the original name does not prove failure: an earlier argument
  # may already have forced this promise. Force the original promise here, at
  # the incumbent boundary, so a cached value remains valid. This handler is
  # restricted to the proved direct, unbound symbol; user factories and bound
  # active/delayed promises never enter it.
  tryCatch({ value; invisible(NULL) }, error = function(e)
    stop(conditionMessage(e), ". ",
         "`data=` is supported only with the formula interface. For native calls, supply vectors, matrices or data frames directly, or use `with(data, ...)`.",
         call. = FALSE))
}
