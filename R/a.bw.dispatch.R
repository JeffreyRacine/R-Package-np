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
