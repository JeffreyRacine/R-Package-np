# Native owners signal this class only when an existing terminal ZERO_RADIUS
# status is being raised. Names/context remain lazy on every successful call.
.np_with_nn_radius_context <- function(expr, continuous.names, context = NULL,
                                       bws = NULL, ntrain = NULL,
                                       leave.one.out = FALSE) {
  withCallingHandlers(
    expr,
    np_nn_zero_radius = function(e) {
      coordinate <- e[["continuous.index", exact = TRUE]]
      if (is.null(e[["variable", exact = TRUE]]) &&
          length(coordinate) == 1L && !is.na(coordinate) &&
          coordinate >= 1L && coordinate <= length(continuous.names)) {
        variable <- continuous.names[[coordinate]]
        if (!is.na(variable) && nzchar(variable)) {
          e[["variable"]] <- variable
          e[["message"]] <- sprintf("continuous variable '%s': %s",
                                      variable, e[["reason", exact = TRUE]])
        }
      }
      if (!is.null(context)) {
        e[["stage"]] <- context
        e[["message"]] <- paste0(context, ": ", conditionMessage(e))
      }
      stop(e)
    },
    error = function(e) {
      # Diagnose only an existing terminal bandwidth/helper failure. Reuse the
      # canonical admission rule, and leave all metadata lazy on success.
      if (!(conditionMessage(e) %in% c(
          "\n** Error: invalid bandwidth.",
          "C_np_regression_lp_apply_conditional: LP hat helper failed",
          "C_np_regression_lp_apply_conditional: LP apply helper failed")) ||
          is.null(bws) || is.null(ntrain))
        return(invisible(NULL))
      effective.n <- ntrain - as.integer(isTRUE(leave.one.out))
      failure <- tryCatch(
        npValidateExtendedNnContinuousBandwidth(bws,
          where = sprintf("effective training size %d", effective.n),
          nobs = effective.n),
        error = identity)
      if (inherits(failure, "error")) {
        e[["message"]] <- paste0(conditionMessage(e), "; ",
                                  conditionMessage(failure))
        stop(e)
      }
    }
  )
}
