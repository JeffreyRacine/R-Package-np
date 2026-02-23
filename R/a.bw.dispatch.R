.np_bw_dispatch_target <- function(dots, data_arg_names = character(), eval_env = parent.frame()) {
  if (length(dots) == 0L)
    stop("invoked without arguments")

  dot.names <- names(dots)

  if (!is.null(dot.names) && any(dot.names == "formula"))
    return(eval(dots[[which(dot.names == "formula")[1L]]], envir = eval_env))

  first.val <- eval(dots[[1L]], envir = eval_env)
  if (inherits(first.val, "formula"))
    return(first.val)

  if (!is.null(dot.names) && any(dot.names == "bws"))
    return(eval(dots[[which(dot.names == "bws")[1L]]], envir = eval_env))

  if (length(data_arg_names) > 0L &&
      !is.null(dot.names) &&
      any(dot.names %in% data_arg_names))
    return(NULL)

  first.val
}

.np_bw_formula_from_call <- function(call_obj, eval_env = parent.frame()) {
  if (missing(call_obj) || !is.call(call_obj))
    return(NULL)

  for (i in seq_along(call_obj)) {
    val <- tryCatch(eval(call_obj[[i]], envir = eval_env), error = function(e) NULL)
    if (inherits(val, "formula"))
      return(call_obj[[i]])
  }

  NULL
}

.np_bw_resolve_formula <- function(formula_obj, formula_call = NULL, eval_env = parent.frame()) {
  if (is.null(formula_call))
    return(formula_obj)

  resolved <- tryCatch(eval(formula_call, envir = eval_env), error = function(e) NULL)
  if (inherits(resolved, "formula"))
    return(resolved)

  formula_obj
}
