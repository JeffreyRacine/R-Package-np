# Validate the effective public outcome without coercing, copying or scanning it.
# Predictor scores and operator right-hand sides have separate internal owners.
.np_require_numeric_response <- function(y, where, argument) {
  if (!is.numeric(y) || !(is.integer(y) || is.double(y)) ||
      is.factor(y) || !is.null(dim(y)))
    .np_numeric_response_error(where, argument)
  invisible(NULL)
}

.np_numeric_response_error <- function(where, argument) {
  stop(paste0(where, "(): '", argument,
              "' must be a numeric (integer or double) response vector; ",
              "factors and ordered factors are not supported."), call. = FALSE)
}

# Reconstruction must not erase a known categorical response on an older fit.
# Absent metadata cannot establish the original type; preserve that legacy path.
.np_plot_reconstruct_response <- function(object, bws, where) {
  type <- bws[["ydati", exact = TRUE]]
  if (is.list(type) &&
      any(c(type[["iuno", exact = TRUE]], type[["iord", exact = TRUE]]),
          na.rm = TRUE))
    .np_numeric_response_error(where, "ydat")
  object[["mean", exact = TRUE]] + object[["resid", exact = TRUE]]
}

# Plot engines resolve public outcomes independently of the estimator entries.
# Check that raw value before NA indexing or factor scoring can change its shape.
.np_plot_require_response <- function(y, bws, where) {
  if (inherits(bws, "sibandwidth") &&
      identical(bws[["method", exact = TRUE]], "kleinspady"))
    .npindex_check_binary_response(y, paste0(where, "() 'ydat'"))
  else
    .np_require_numeric_response(y, where, "ydat")
  invisible(NULL)
}

# Forward only supplied arguments through their existing promises. The private
# context records call semantics; it is neither an object flag nor a public option.
.np_regression_forward <- function(implementation, frame, call, expanded,
                                   caller, response, dots, public) {
  matched <- match.call(definition = implementation, call = expanded,
                        expand.dots = FALSE)
  context <- list(call = matched, expanded = expanded, caller = caller,
                  response = response, dots = dots, public = public)
  args <- as.list(call)[-1L]
  forwarded <- list()
  for (name in names(args)) {
    if (identical(name, "...")) {
      values <- as.list(args[[name]])
      for (i in seq_along(values)) values[[i]] <- as.name(paste0("..", i))
      forwarded <- c(forwarded, values)
    } else forwarded[[name]] <- as.name(name)
  }
  entry <- new.env(parent = frame)
  entry$.np.context <- context
  forwarded <- c(forwarded, list(.np.context = quote(.np.context)))
  invocation <- as.call(c(list(implementation), forwarded))
  public.call <- sys.call(sys.parent())
  withCallingHandlers(eval(invocation, envir = entry), error = function(condition) {
    # Errors raised by the implementation belong to the original entry call.
    # Helper-owned conditions, classes, messages and null-call errors propagate
    # unchanged. Re-signal once; never recover or continue a failed computation.
    if (identical(conditionCall(condition), invocation)) {
      condition$call <- public.call
      stop(condition)
    }
  })
}

# This entry is used only for the response-on-Z children owned by npplreg.
# The existing selector and numerical implementations remain the owners.
.npregbw_auxiliary <- function(...) {
  call <- match.call(expand.dots = FALSE)
  target <- .np_bw_dispatch_target(dots = call$...,
    data_arg_names = c("xdat", "ydat"), eval_env = parent.frame(),
    promise.frame = environment())
  implementation <- if (is.null(target)) .npregbw_null_impl else
    if (inherits(target, "rbandwidth")) .npregbw_bandwidth_impl else
      .npregbw_default_impl
  expanded <- match.call()
  call[[1L]] <- expanded[[1L]] <- if (is.null(target)) quote(npregbw.NULL) else
    if (inherits(target, "rbandwidth")) quote(npregbw.rbandwidth) else
      quote(npregbw.default)
  matched <- match.call(definition = implementation, call = expanded,
                        expand.dots = FALSE)
  .np_regression_forward(implementation, environment(), call, expanded,
    parent.frame(), if (is.null(matched$ydat)) quote(ydat) else matched$ydat,
    substitute(list(...))[-1L], FALSE)
}
