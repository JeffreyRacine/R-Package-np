# Missing uncertainty is an extraction condition, never an implicit refit.
.np_se_refit_hint <- function(expr, family, switches = "se = TRUE",
                              bws.field = "bws", data.hint = NULL) {
  object <- if (is.symbol(expr) && nzchar(as.character(expr))) {
    paste(deparse(expr, backtick = TRUE), collapse = "")
  } else {
    "object"
  }
  paste0("Refit without repeating bandwidth search: ", family, "(bws = ",
         object, "$", bws.field, ", ", switches,
         if (!is.null(data.hint)) paste0(", ", data.hint) else "", ").",
         if (identical(object, "object"))
           " Replace 'object' with your fitted model." else "")
}

.np_stop_missing_output <- function(expr, family,
                                    message = "standard errors were not computed.",
                                    switches = "se = TRUE",
                                    bws.field = "bws", data.hint = NULL) {
  stop(paste(message,
             .np_se_refit_hint(expr, family, switches, bws.field, data.hint)),
       call. = FALSE)
}

# Numeric NA is a computed undefined row, not the legacy logical NA placeholder.
.np_se_output_missing <- function(value) {
  is.null(value) || length(value) == 0L ||
    (is.logical(value) && length(value) == 1L && is.na(value))
}

.np_require_stored_se <- function(x, value, family, what = "standard errors",
                                   expr = substitute(x), switches = "se = TRUE",
                                   bws.field = "bws", data.hint = NULL) {
  if (identical(x[["se", exact = TRUE]], FALSE) ||
      .np_se_output_missing(value))
    .np_stop_missing_output(expr, family, paste0(what, " were not computed."),
                            switches, bws.field, data.hint)
  value
}
