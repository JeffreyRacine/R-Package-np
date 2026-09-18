# Native bandwidth-object fits share the evaluation roles already supported by
# fitted-object predict methods. Do not consult a call or an ambient frame.
.np_native_newdata_parts <- function(newdata, groups, where) {
  nd <- toFrame(newdata)
  if (length(groups) == 1L)
    return(setNames(list(nd), names(groups)))
  required <- unlist(groups, use.names = FALSE)
  if (any(lengths(groups) == 0L) || anyNA(required) || any(!nzchar(required)) ||
      anyDuplicated(required) || anyDuplicated(names(nd)))
    stop(sprintf("%s: native 'newdata' roles are ambiguous; supply %s explicitly",
                 where, paste(shQuote(names(groups)), collapse = " and ")),
         call. = FALSE)
  if (!all(required %in% names(nd)))
    stop(sprintf("%s: 'newdata' must include columns %s, or supply %s explicitly",
                 where, paste(shQuote(required), collapse = ", "),
                 paste(shQuote(names(groups)), collapse = " and ")),
         call. = FALSE)
  lapply(groups, function(columns) nd[, columns, drop = FALSE])
}
