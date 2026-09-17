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
    get_dot_value <- function(env) {
      dots <- tryCatch(
        eval(substitute(list(...)), envir = env),
        error = function(e) not_found
      )
      if (identical(dots, not_found) || !is.list(dots))
        return(not_found)
      dot_names <- names(dots)
      if (is.null(dot_names) || !(sym %in% dot_names))
        return(not_found)
      dots[[which(dot_names == sym)[1L]]]
    }
    eval_symbol <- function(env) {
      tryCatch(eval(expr, envir = env), error = function(e) not_found)
    }

    sym_val <- get0(sym, envir = eval_env, inherits = TRUE, ifnotfound = not_found)
    if (!identical(sym_val, not_found))
      return(list(ok = TRUE, value = sym_val, error = NULL))
    sym_val <- get_dot_value(eval_env)
    if (!identical(sym_val, not_found))
      return(list(ok = TRUE, value = sym_val, error = NULL))
    sym_val <- eval_symbol(eval_env)
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
        sym_val <- get_dot_value(env_i)
        if (!identical(sym_val, not_found))
          return(list(ok = TRUE, value = sym_val, error = NULL))
        sym_val <- eval_symbol(env_i)
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
    function(env) tryCatch(
      list(ok = TRUE, value = eval(expr, envir = env), error = NULL),
      error = function(e) list(ok = FALSE, value = NULL, error = e))
  } else {
    # Same contract path when an explicit enclos is supplied.
    function(env) tryCatch(
      list(ok = TRUE, value = eval(expr, envir = env, enclos = enclos), error = NULL),
      error = function(e) list(ok = FALSE, value = NULL, error = e))
  }

  val <- eval_once(eval_env)
  if (isTRUE(val$ok))
    return(val)

  first_error <- val$error
  if (!isTRUE(search_frames))
    return(list(ok = FALSE, value = NULL, error = first_error))

  frames <- sys.frames()
  for (i in rev(seq_along(frames))) {
    env_i <- frames[[i]]
    if (identical(env_i, eval_env))
      next
    val_i <- eval_once(env_i)
    if (isTRUE(val_i$ok))
      return(val_i)
  }

  list(ok = FALSE, value = NULL, error = first_error)
}

.np_formula_expand_call <- function(call_obj, caller_env = parent.frame()) {
  if (any(vapply(as.list(call_obj)[-1L], identical, logical(1L), quote(...)))) {
    call_obj <- match.call(definition = function(...) NULL, call = call_obj,
                           expand.dots = TRUE, envir = caller_env)
    # Keep ordinary arguments as references to the original dot promises.
    # subset alone is syntax evaluated in the model-frame data mask; ..n is
    # not meaningful in that mask or in a separately constructed formula env.
    subset <- match("subset", names(call_obj), nomatch = 0L)
    if (subset) call_obj[[subset]] <-
      .np_formula_dot_expression(call_obj[[subset]], caller_env)
  }
  call_obj
}

.np_formula_dot_expression <- function(expr, caller_env) {
  if (is.symbol(expr) && grepl("^\\.\\.[0-9]+$", as.character(expr))) {
    index <- as.integer(substring(as.character(expr), 3L))
    expressions <- eval(quote(substitute(list(...))), envir = caller_env)
    if (index > 0L && index < length(expressions)) return(expressions[[index + 1L]])
  }
  expr
}

.np_formula_value <- function(formula = NULL, bws, data, data.name = "xdat") {
  if (inherits(formula, "formula"))
    return(list(formal = "formula", value = formula, data.name = data.name))
  if (!missing(bws) && inherits(bws, "formula"))
    return(list(formal = "bws", value = bws, data.name = data.name))
  if (!missing(data) && inherits(data, "formula"))
    return(list(formal = data.name, value = data, data.name = data.name))
  NULL
}

.np_eval_bw_call <- function(call_obj, caller_env = parent.frame(),
                             formula.value = NULL) {
  if (!is.call(call_obj))
    stop("bandwidth selector call is malformed", call. = FALSE)

  # Namespace-only estimator calls reconstruct an unqualified selector name.
  # Resolve that name before executing; never replay a failed search in other
  # frames. The existence check must not force a caller's active binding, and
  # arguments must still be evaluated in the original caller environment.
  selector <- call_obj[[1L]]
  if (is.symbol(selector) &&
      !exists(as.character(selector), envir = caller_env, inherits = TRUE) &&
      as.character(selector) %in% getNamespaceExports("np"))
    call_obj[[1L]] <- call("::", as.name("np"), selector)

  formula.expression <- NULL
  if (!is.null(formula.value)) {
    if (!inherits(formula.value$value, "formula"))
      stop("invalid resolved formula handoff", call. = FALSE)
    # Match positions without evaluating expressions. The selector handoff
    # retains the fitter's bws/data positions, including a positional formula
    # beside named numeric bws. Do not guess from the first unnamed argument:
    # an owner can intentionally discard its formula for explicit native data.
    tagged <- call_obj
    for (i in seq.int(2L, length(tagged))) tagged[[i]] <- i
    matcher <- function(bws, xdat, ...) NULL
    names(formals(matcher))[2L] <- formula.value$data.name
    matched <- match.call(matcher, tagged, expand.dots = TRUE)
    source <- formula.value$formal
    if (source == "bws" && "formula" %in% names(matched)) source <- "formula"
    index <- matched[[source]]
    if (!is.null(index)) {
      formula.expression <- .np_formula_dot_expression(call_obj[[index]], caller_env)
      call_obj[[index]] <- substitute(quote(VALUE), list(VALUE = formula.value$value))
    }
  }
  result <- eval(call_obj, envir = caller_env)
  # Execute with the already-resolved value, but preserve the user's formula
  # expression in the existing call metadata. No handoff state is retained.
  if (!is.null(formula.expression) && !is.null(result[["formula"]]) &&
      is.call(result[["call"]]))
    result[["call"]][["formula"]] <- formula.expression
  result
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
  if (!is.null(nomad)) {
    nomad.mode <- tryCatch(
      npValidateNomadControl(nomad, "nomad"),
      error = function(e) "false"
    )
    if (nomad.mode %in% c("true", "auto"))
      return(TRUE)
  }

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

.np_bw_dispatch_target <- function(dots, data_arg_names = character(), eval_env = parent.frame(),
                                   promise.frame = NULL) {
  if (length(dots) == 0L)
    stop("invoked without arguments")

  dot.names <- names(dots)
  has.named.bws <- !is.null(dot.names) && any(dot.names == "bws")

  value <- function(i) {
    if (is.null(promise.frame))
      return(.np_try_eval_in_frames(dots[[i]], eval_env = eval_env))
    # Force the original dot promise, not its captured syntax. UseMethod
    # forwards this same promise to the selected method. In particular an
    # error must propagate here, not cause the formula factory to be retried.
    list(ok = TRUE, value = eval(substitute(...elt(INDEX), list(INDEX = i)),
                                envir = promise.frame))
  }

  if (!is.null(dot.names) && any(dot.names == "formula")) {
    fval <- value(which(dot.names == "formula")[1L])
    if (isTRUE(fval$ok))
      return(fval$value)
  }

  first.eval <- value(1L)
  if (!isTRUE(first.eval$ok))
    return(NULL)
  first.val <- first.eval$value
  if (inherits(first.val, "formula"))
    return(first.val)

  if (!has.named.bws)
    return(NULL)

  if (has.named.bws) {
    bval <- value(which(dot.names == "bws")[1L])
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

# Used inside model.frame's predvars: formula expressions are evaluated once.
# Ordinary columns already belong to the intersection; model.frame checks their
# lengths and applies subset/na.action after the time-series alignment.
.np_formula_align_values <- function(values) {
  is.ts <- vapply(values, inherits, logical(1), "ts")
  if (!any(is.ts)) return(values)
  series <- values[is.ts]
  if (length(series) > 1L) {
    widths <- vapply(series, NCOL, integer(1))
    # ts.intersect expands mts columns. Reassemble the original expressions
    # before model.frame assigns its variable names; never recycle columns
    # back into expression slots.
    inputs <- series
    names(inputs) <- paste0(".series", seq_along(inputs))
    aligned <- do.call(stats::ts.intersect, inputs)
    if (is.null(aligned))
      stop("time-series formula variables have no common observations", call. = FALSE)
    if (NCOL(aligned) != sum(widths))
      stop("time-series formula alignment changed the number of columns", call. = FALSE)
    offsets <- c(0L, cumsum(widths))
    series <- lapply(seq_along(series), function(i) {
      columns <- offsets[[i]] + seq_len(widths[[i]])
      if (is.matrix(series[[i]])) {
        value <- aligned[, columns, drop = FALSE]
        colnames(value) <- colnames(series[[i]])
      } else {
        value <- aligned[, columns]
      }
      value
    })
  }
  values[is.ts] <- lapply(series, function(x) {
    attr(x, "tsp") <- NULL
    class(x) <- setdiff(class(x), c("ts", "mts"))
    x
  })
  values
}

.np_formula_aligned_terms <- function(tt) {
  .np_formula_validate_terms(tt)
  # A namespace-qualified prediction expression survives saved-bandwidth
  # refits/newdata without capturing a caller frame or the training data.
  attr(tt, "predvars") <- substitute(
    utils::getFromNamespace(".np_formula_align_values", PACKAGE)(VARIABLES),
    list(PACKAGE = "np", VARIABLES = attr(tt, "variables")))
  tt
}

# Conditional bandwidths retain expanded, one-sided joint terms, not ordinary
# response-bearing terms. Drop their response before rebuilding prediction
# expressions; reparsing the original formula would lose a fitted `.` expansion.
.np_formula_conditional_rhs_terms <- function(bws) {
  .np_formula_validate_syntax(bws$formula, conditional.response = TRUE)
  tt <- bws$terms
  attr(tt, "predvars") <- .np_formula_unwrap_prediction(tt)$prediction
  response <- match(bws$variableNames[["response"]], attr(tt, "term.labels"))
  drop.terms(tt, response)
}
