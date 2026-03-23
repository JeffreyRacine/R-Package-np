.np_missing_binding_sentinel <- new.env(parent = emptyenv())

.np_try_eval_in_frames <- function(expr, eval_env = parent.frame(), enclos = NULL, search_frames = TRUE) {
  if (!is.language(expr))
    return(list(ok = TRUE, value = expr, error = NULL))

  sym <- NULL
  not_found <- .np_missing_binding_sentinel
  if (is.symbol(expr)) {
    sym <- as.character(expr)
  }

  if (!is.null(sym) && is.environment(eval_env)) {
    sym_val <- get0(sym, envir = eval_env, inherits = TRUE, ifnotfound = not_found)
    if (!identical(sym_val, not_found))
      return(list(ok = TRUE, value = sym_val, error = NULL))

    first_error <- simpleError(sprintf("object '%s' not found", sym))
    if (!isTRUE(search_frames))
      return(list(ok = FALSE, value = NULL, error = first_error))

    frames <- sys.frames()
    for (i in rev(seq_along(frames))) {
      env_i <- frames[[i]]
      if (identical(env_i, eval_env))
        next
      if (is.environment(env_i)) {
        sym_val <- get0(sym, envir = env_i, inherits = TRUE, ifnotfound = not_found)
        if (!identical(sym_val, not_found))
          return(list(ok = TRUE, value = sym_val, error = NULL))
      }
    }

    return(list(ok = FALSE, value = NULL, error = first_error))
  }

  eval_once <- if (is.null(enclos)) {
    # Intentional dynamic evaluation: selector/formula NSE semantics require
    # frame-aware call evaluation across caller frames.
    # Guarded by helper contracts in tests/testthat:
    # - test-call-eval-helpers-contract.R
    # - test-bw-eval-helper-contract.R
    function(env) tryCatch(eval(expr, envir = env), error = function(e) e)
  } else {
    # Same contract path when an explicit enclos is supplied.
    function(env) tryCatch(eval(expr, envir = env, enclos = enclos), error = function(e) e)
  }

  val <- eval_once(eval_env)
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
    val_i <- eval_once(env_i)
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

.np_bw_call_uses_nomad_degree_search <- function(call_obj, caller_env = parent.frame()) {
  if (!is.call(call_obj) || length(call_obj) < 1L)
    return(FALSE)

  fn_expr <- call_obj[[1L]]
  fn <- .np_try_eval_in_frames(fn_expr, eval_env = caller_env)
  if (!isTRUE(fn$ok) || !is.function(fn$value))
    return(FALSE)

  fn_def <- fn$value
  fn_name <- if (is.symbol(fn_expr)) as.character(fn_expr) else NULL
  if (!is.null(fn_name) && nzchar(fn_name)) {
    default_name <- paste0(fn_name, ".default")
    default_fn <- get0(default_name, envir = environment(fn_def), inherits = TRUE)
    if (is.function(default_fn))
      fn_def <- default_fn
  }

  matched <- tryCatch(
    match.call(definition = fn_def, call = call_obj, expand.dots = FALSE),
    error = function(e) NULL
  )
  if (is.null(matched))
    return(FALSE)

  defaults <- formals(fn_def)
  is_missing_arg <- function(z) missing(z)
  arg_expr <- function(source, name) {
    if (is.null(source) || is.null(names(source)) || !(name %in% names(source)))
      return(NULL)
    expr <- source[[name]]
    if (is_missing_arg(expr))
      return(NULL)
    expr
  }
  arg_value <- function(name) {
    expr <- arg_expr(matched, name)
    if (is.null(expr))
      expr <- arg_expr(defaults, name)
    if (is.null(expr))
      return(NULL)
    value <- .np_try_eval_in_frames(expr, eval_env = caller_env)
    if (!isTRUE(value$ok))
      return(NULL)
    value$value
  }

  nomad <- arg_value("nomad")
  if (!is.null(nomad) && isTRUE(as.logical(nomad)[1L]))
    return(TRUE)

  regtype <- arg_value("regtype")
  if (is.null(regtype) || !identical(as.character(regtype)[1L], "lp"))
    return(FALSE)

  degree.select <- arg_value("degree.select")
  if (is.null(degree.select))
    return(FALSE)
  degree.select <- as.character(degree.select)[1L]
  if (!nzchar(degree.select) || identical(degree.select, "manual"))
    return(FALSE)

  search.engine <- arg_value("search.engine")
  if (is.null(search.engine))
    return(FALSE)
  search.engine <- as.character(search.engine)[1L]

  search.engine %in% c("nomad", "nomad+powell")
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
