# Native bandwidth-object fits share the evaluation roles already supported by
# fitted-object predict methods. Do not consult a call or an ambient frame.
.np_args_with_defaults <- function(defaults, explicit) {
  c(defaults[!names(defaults) %in% names(explicit)], explicit)
}

# Match positions/partial names before retained training data fills omissions.
# Match integer tokens, not argument values: language objects remain data and
# no supplied expression is evaluated again. Prediction reserves training slots.
.np_match_native_args <- function(args, definition, retained = character()) {
  if (!length(args)) return(args)
  nms <- names(args)
  formal.names <- names(formals(definition))
  if (!is.null(nms) && all(nzchar(nms)) && all(nms %in% formal.names))
    return(args)
  tokens <- as.list(seq_along(args))
  names(tokens) <- nms
  match.tokens <- function(x) as.list(match.call(definition,
    as.call(c(list(as.name("native"), bws = 0L), x)),
    expand.dots = TRUE))[-1L]
  if (length(retained)) {
    named <- if (is.null(nms)) rep(FALSE, length(args)) else nzchar(nms)
    supplied <- names(match.tokens(tokens[named]))
    reserve <- setdiff(retained, supplied)
    tokens <- c(setNames(rep(list(0L), length(reserve)), reserve), tokens)
  }
  matched <- match.tokens(tokens)
  positions <- unlist(matched, use.names = FALSE)
  keep <- positions > 0L
  setNames(args[positions[keep]], names(matched)[keep])
}

.np_retained_training_args <- function(bws, roles, explicit, definition,
                                       training = names(roles)) {
  explicit <- .np_match_native_args(explicit, definition)
  # Response-only replacement retains the design (e.g. a wild bootstrap).
  # A replaced design must still supply every required training role.
  if (!any(setdiff(training, "tydat") %in% names(explicit))) {
    for (arg in setdiff(names(roles), names(explicit)))
      explicit[arg] <- list(.np_eval_bws_call_arg(bws, roles[[arg]]))
  } else {
    # A missing dispatch role would otherwise re-enter this call method.
    # Preserve concrete-method required roles and optional defaults (e.g. z=x).
    for (arg in setdiff(training, names(explicit))) {
      default <- formals(definition)[[arg]]
      if (is.call(default) && identical(default[[1L]], quote(stop)))
        stop(sprintf("training data '%s' missing", arg), call. = FALSE)
    }
  }
  .np_args_with_defaults(list(bws = bws), explicit)
}

.np_native_newdata_parts <- function(newdata, groups, where) {
  # Preserve the established positional input for an unnamed single-role
  # matrix/vector. Once names are supplied, they identify variables, not rows.
  named <- is.data.frame(newdata) ||
    (!is.null(dim(newdata)) && !is.null(colnames(newdata)))
  nd <- toFrame(newdata)
  if (length(groups) == 1L && !named)
    return(setNames(list(nd), names(groups)))
  required <- unlist(groups, use.names = FALSE)
  if (any(lengths(groups) == 0L) || anyNA(required) || any(!nzchar(required)) ||
      anyDuplicated(required) || is.null(names(nd)) || anyNA(names(nd)) ||
      any(!nzchar(names(nd))) || anyDuplicated(names(nd)))
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
