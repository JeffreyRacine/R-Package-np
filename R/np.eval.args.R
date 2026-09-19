# Native bandwidth-object fits share the evaluation roles already supported by
# fitted-object predict methods. Do not consult a call or an ambient frame.
.np_args_with_defaults <- function(defaults, explicit) {
  c(defaults[!names(defaults) %in% names(explicit)], explicit)
}

.np_retained_training_args <- function(bws, roles, explicit) {
  for (arg in setdiff(names(roles), names(explicit)))
    explicit[arg] <- list(.np_eval_bws_call_arg(bws, roles[[arg]]))
  .np_args_with_defaults(list(bws = bws), explicit)
}

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
