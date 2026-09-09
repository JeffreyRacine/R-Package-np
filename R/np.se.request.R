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

.np_require_stored_se <- function(x, value, family, what = "standard errors",
                                   expr = substitute(x), switches = "se = TRUE",
                                   bws.field = "bws", data.hint = NULL) {
  if (identical(x[["se", exact = TRUE]], FALSE) ||
      is.null(value) || length(value) == 0L)
    stop(paste0(what, " were not computed. ",
                .np_se_refit_hint(expr, family, switches, bws.field, data.hint)),
         call. = FALSE)
  value
}
