.np_try_eval_in_frames <- function(expr, eval_env = parent.frame(), enclos = NULL, search_frames = TRUE) {
  if (is.symbol(expr) &&
      is.environment(eval_env) &&
      exists(as.character(expr), envir = eval_env, inherits = TRUE)) {
    return(list(ok = TRUE, value = get(as.character(expr), envir = eval_env, inherits = TRUE), error = NULL))
  }

  val <- tryCatch(
    if (is.null(enclos)) eval(expr, envir = eval_env) else eval(expr, envir = eval_env, enclos = enclos),
    error = function(e) e
  )
  if (!inherits(val, "error"))
    return(list(ok = TRUE, value = val, error = NULL))

  first_error <- val
  if (!isTRUE(search_frames))
    return(list(ok = FALSE, value = NULL, error = first_error))

  frames <- sys.frames()
  for (i in rev(seq_along(frames))) {
    env_i <- frames[[i]]
    if (identical(env_i, eval_env))
      next
    if (is.symbol(expr) &&
        is.environment(env_i) &&
        exists(as.character(expr), envir = env_i, inherits = TRUE)) {
      return(list(ok = TRUE, value = get(as.character(expr), envir = env_i, inherits = TRUE), error = NULL))
    }
    val_i <- tryCatch(
      if (is.null(enclos)) eval(expr, envir = env_i) else eval(expr, envir = env_i, enclos = enclos),
      error = function(e) e
    )
    if (!inherits(val_i, "error"))
      return(list(ok = TRUE, value = val_i, error = NULL))
  }

  list(ok = FALSE, value = NULL, error = first_error)
}

.np_eval_bw_call <- function(call_obj, caller_env = parent.frame()) {
  if (!is.call(call_obj))
    stop("bandwidth selector call is malformed", call. = FALSE)

  fallback <- .np_try_eval_in_frames(call_obj, eval_env = caller_env)
  if (isTRUE(fallback$ok))
    return(fallback$value)

  if (inherits(fallback$error, "error"))
    stop(conditionMessage(fallback$error), call. = FALSE)
  stop("bandwidth selector call evaluation failed", call. = FALSE)
}

.np_bw_dispatch_target <- function(dots, data_arg_names = character(), eval_env = parent.frame()) {
  if (length(dots) == 0L)
    stop("invoked without arguments")

  dot.names <- names(dots)
  has.named.bws <- !is.null(dot.names) && any(dot.names == "bws")
  has.named.data <- length(data_arg_names) > 0L &&
    !is.null(dot.names) &&
    any(dot.names %in% data_arg_names)

  if (!is.null(dot.names) && any(dot.names == "formula")) {
    fval <- .np_try_eval_in_frames(dots[[which(dot.names == "formula")[1L]]], eval_env = eval_env)
    if (isTRUE(fval$ok))
      return(fval$value)
  }

  if (has.named.data && !has.named.bws)
    return(NULL)

  first.eval <- .np_try_eval_in_frames(dots[[1L]], eval_env = eval_env)
  if (!isTRUE(first.eval$ok))
    return(NULL)
  first.val <- first.eval$value
  if (inherits(first.val, "formula"))
    return(first.val)

  if (has.named.bws) {
    bval <- .np_try_eval_in_frames(dots[[which(dot.names == "bws")[1L]]], eval_env = eval_env)
    if (isTRUE(bval$ok))
      return(bval$value)
  }

  first.val
}

.np_bw_formula_from_call <- function(call_obj, eval_env = parent.frame()) {
  if (missing(call_obj) || !is.call(call_obj))
    return(NULL)

  if (length(call_obj) < 2L)
    return(NULL)

  for (i in 2:length(call_obj)) {
    val <- .np_try_eval_in_frames(call_obj[[i]], eval_env = eval_env)
    if (isTRUE(val$ok) && inherits(val$value, "formula"))
      return(call_obj[[i]])
  }

  NULL
}

.np_bw_resolve_formula <- function(formula_obj, formula_call = NULL, eval_env = parent.frame()) {
  if (is.null(formula_call))
    return(formula_obj)

  resolved <- .np_try_eval_in_frames(formula_call, eval_env = eval_env)
  if (isTRUE(resolved$ok) && inherits(resolved$value, "formula"))
    return(resolved$value)

  formula_obj
}

.np_terms_variable_values <- function(terms_obj, data, eval_env = environment(terms_obj)) {
  if (missing(terms_obj) || is.null(terms_obj))
    return(list())

  vars <- attr(terms_obj, "variables")
  if (is.null(vars))
    return(list())

  out <- .np_try_eval_in_frames(vars, eval_env = data, enclos = eval_env)
  if (!isTRUE(out$ok) || is.null(out$value))
    return(list())
  out <- out$value
  if (!is.list(out))
    out <- as.list(out)

  out
}

.np_terms_ts_mask <- function(terms_obj, data, eval_env = environment(terms_obj)) {
  vals <- .np_terms_variable_values(terms_obj = terms_obj, data = data, eval_env = eval_env)
  if (length(vals) == 0L)
    return(logical(0))

  vapply(vals, inherits, logical(1), "ts")
}
