# Native owners signal this class only when an existing terminal ZERO_RADIUS
# status is being raised. Names/context remain lazy on every successful call.
.np_with_nn_radius_context <- function(expr, continuous.names, context = NULL) {
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
    }
  )
}
