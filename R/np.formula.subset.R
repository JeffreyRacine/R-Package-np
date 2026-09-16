# Evaluate expression values once, align them transiently, and retain portable
# prediction expressions (including makepredictcall metadata), never the data.
.np_formula_unwrap_prediction <- function(tt) {
  prediction <- attr(tt, "predvars")
  variables <- attr(tt, "variables")
  owned <- is.call(prediction) && length(prediction) == 2L &&
    is.call(prediction[[1L]]) && length(prediction[[1L]]) == 3L &&
    identical(prediction[[1L]][[1L]], quote(utils::getFromNamespace)) &&
    identical(prediction[[1L]][[2L]], ".np_formula_align_values") &&
    is.character(prediction[[1L]][[3L]]) &&
    prediction[[1L]][[3L]] %in% c("np", "npRmpi")
  if (owned) prediction <- prediction[[2L]]
  # Older conditional objects used base-R alignment scaffolding. Recognize
  # only exact reconstructions from the stored variable expressions: arbitrary
  # user prediction calls and their trained metadata are not ours to replace.
  legacy <- FALSE
  if (is.call(prediction) && identical(prediction[[1L]], quote(as.data.frame))) {
    expected <- as.call(list(quote(as.data.frame),
      as.call(c(list(quote(ts.intersect)), as.list(variables)[-1L]))))
    legacy <- identical(prediction, expected)
  } else if (is.call(prediction) && identical(prediction[[1L]], quote(`[`)) &&
             length(prediction) == 4L && is.call(prediction[[2L]]) &&
             identical(prediction[[2L]][[1L]], quote(cbind)) &&
             length(prediction[[2L]]) >= 2L) {
    frame <- prediction[[2L]][[2L]]
    if (is.call(frame) && identical(frame[[1L]], quote(as.data.frame)) &&
        length(frame) == 2L && is.call(frame[[2L]]) &&
        identical(frame[[2L]][[1L]], quote(ts.intersect))) {
      arguments <- as.list(variables)[-1L]
      series <- as.list(frame[[2L]])[-1L]
      is.ts <- vapply(arguments, function(x)
        any(vapply(series, identical, logical(1L), x)), logical(1L))
      if (any(is.ts) && any(!is.ts) && sum(is.ts) == length(series)) {
        combined <- as.call(c(list(quote(cbind),
          as.call(list(quote(as.data.frame),
            as.call(c(list(quote(ts.intersect)), arguments[is.ts]))))),
          arguments[!is.ts], list(check.rows = TRUE)))
        index <- order(c(which(is.ts), which(!is.ts)))
        expected <- substitute(COMBINED[, INDEX], list(COMBINED = combined, INDEX = index))
        legacy <- identical(prediction, expected)
      }
    }
  }
  if (legacy) prediction <- variables
  list(prediction = if (is.null(prediction)) attr(tt, "variables") else prediction,
       makepredictcall = is.null(attr(tt, "predvars")) || owned || legacy)
}

.np_formula_model_frame <- function(formula, data = NULL, subset, na.action,
                                    drop.unused.levels = FALSE, xlev = NULL, ...) {
  if (!is.data.frame(data) && !is.environment(data) && !is.null(attr(data, "class")))
    data <- as.data.frame(data)
  tt <- if (inherits(formula, "terms")) formula else terms(formula, data = data)
  prepared <- .np_formula_unwrap_prediction(tt)
  prediction <- prepared$prediction
  values <- eval(prediction, data, environment(tt))
  if (prepared$makepredictcall) {
    variables <- attr(tt, "variables")
    for (i in seq_along(values))
      prediction[[i + 1L]] <- stats::makepredictcall(values[[i]], variables[[i + 1L]])
  }
  attr(tt, "predvars") <- .np_formula_align_values(values)
  frame.call <- match.call()
  frame.call[[1L]] <- quote(stats::model.frame)
  frame.call[["formula"]] <- tt
  frame.call["data"] <- list(data)
  frame <- eval(frame.call, parent.frame())
  retained <- attr(frame, "terms")
  attr(retained, "predvars") <- prediction
  attr(frame, "terms") <- retained
  frame
}

# An automatic constructor/fit transaction may hand off its frame exactly once.
# The context never travels to native/MPI leaves or remains in a stored call.
.np_formula_frame_store <- function(state, frame) {
  if (is.null(state)) return(invisible(NULL))
  if (!is.environment(state) || !identical(parent.env(state), emptyenv()) ||
      exists("frame", envir = state, inherits = FALSE))
    stop("invalid formula preparation context", call. = FALSE)
  state$frame <- frame
  invisible(NULL)
}

.np_formula_frame_take <- function(state) {
  if (!is.environment(state) || !exists("frame", envir = state, inherits = FALSE))
    stop("formula preparation context has no training frame", call. = FALSE)
  frame <- state$frame
  rm("frame", envir = state)
  frame
}

.np_formula_call_public <- function(mc) {
  mc$.np.formula.state <- NULL
  dots <- mc[["...", exact = TRUE]]
  if (!is.null(dots) && !is.null(names(dots))) {
    keep <- is.na(names(dots)) | names(dots) != ".np.formula.state"
    mc[["..."]] <- dots[keep]
  }
  mc
}

# Resolve formula subsets before any transport or model-frame replay. Return
# values only: the caller environment is used transiently, never retained.
.np_formula_subset_inputs <- function(data, subset.expr, caller) {
  list(data = data, subset = eval(subset.expr, envir = data, enclos = caller))
}

.np_bws_formula_model_frame <- function(bws, mf.args, data.override = FALSE) {
  call.env <- environment(bws$call)
  if (!data.override && is.environment(call.env) &&
      "data" %in% names(mf.args) && is.language(mf.args[["data"]])) {
    # Resolve only the saved data expression in its existing owner, when
    # model.frame forces it. Formula variables and subset retain their lexical
    # data-mask semantics; neither the stored call nor its environment changes.
    mf.args[["data"]] <- substitute(base::eval(quote(EXPR), envir = OWNER),
                                    list(EXPR = mf.args[["data"]],
                                         OWNER = call.env))
  }
  do.call(.np_formula_model_frame, mf.args,
          envir = environment(mf.args[["formula"]]))
}
