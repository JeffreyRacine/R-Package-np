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
  validate <- state[["validate.training"]]
  if (!is.null(validate)) {
    if (!is.function(validate))
      stop("invalid formula training validator", call. = FALSE)
    rm("validate.training", envir = state)
    validate(frame)
  }
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

# Partially linear formulas have two retained role terms but one observation
# set. Keep the joint frame invocation-local and project trained metadata back
# into the existing terms/xterms fields.
.np_plreg_formula_spec <- function(formula, chromoly = NULL) {
  if (is.null(chromoly))
    chromoly <- explodePipe(formula, env = environment(formula))
  if (length(chromoly) != 3L)
    stop("invoked with improper formula, please see npplregbw documentation for proper use")
  bronze <- vapply(chromoly, paste, character(1L), collapse = " + ")
  make <- function(text) terms(as.formula(text, env = environment(formula)))
  list(chromoly = chromoly,
       terms = make(paste(bronze[[1L]], "~", bronze[[3L]])),
       xterms = make(paste("~", bronze[[2L]])),
       joint = make(paste(bronze[[1L]], "~", bronze[[2L]], "+", bronze[[3L]])))
}

.np_plreg_formula_indices <- function(role, joint) {
  variables <- as.list(attr(role, "variables"))[-1L]
  combined <- as.list(attr(joint, "variables"))[-1L]
  vapply(variables, function(variable) {
    index <- which(vapply(combined, identical, logical(1L), variable))
    if (length(index) != 1L)
      stop("inconsistent partially linear formula variables", call. = FALSE)
    index
  }, integer(1L))
}

.np_plreg_formula_prediction <- function(role, joint) {
  prediction <- .np_formula_unwrap_prediction(role)$prediction
  variables <- attr(role, "variables")
  if (is.call(prediction) && identical(prediction[[1L]], quote(list)) &&
      length(prediction) == length(variables))
    return(list(prediction = prediction, scaffold = NULL))

  # Old plreg objects wrapped the full joint alignment in an exact role
  # subset. Recognize that package-owned shape only; user prediction calls
  # and trained transform expressions are never stripped heuristically.
  if (is.call(prediction) && length(prediction) == 5L &&
      identical(prediction[[1L]], quote(`[`))) {
    inner <- prediction[[2L]]
    indices <- .np_plreg_formula_indices(role, joint)
    expected <- substitute(INNER[, INDEX, drop = FALSE],
                           list(INNER = inner, INDEX = indices))
    if (identical(prediction, expected)) {
      aligned <- inner
      if (is.call(aligned) && length(aligned) == 2L &&
          identical(aligned[[1L]], quote(`(`)))
        aligned <- aligned[[2L]]
      temporary <- joint
      attr(temporary, "predvars") <- aligned
      unwrapped <- .np_formula_unwrap_prediction(temporary)$prediction
      if (!identical(aligned, unwrapped) &&
          identical(unwrapped, attr(joint, "variables")))
        return(list(prediction = variables, scaffold = inner))
    }
  }
  stop("unsupported partially linear prediction metadata", call. = FALSE)
}

.np_plreg_formula_terms <- function(bws) {
  spec <- .np_plreg_formula_spec(bws$formula, bws$chromoly)
  joint <- spec$joint
  roles <- list(bws$terms, bws$xterms)
  predictions <- lapply(roles, .np_plreg_formula_prediction, joint = joint)
  if (!is.null(predictions[[1L]]$scaffold) &&
      !is.null(predictions[[2L]]$scaffold) &&
      !identical(predictions[[1L]]$scaffold, predictions[[2L]]$scaffold))
    stop("inconsistent partially linear alignment metadata", call. = FALSE)
  combined <- attr(joint, "variables")
  seen <- rep_len(FALSE, length(combined) - 1L)
  for (i in seq_along(roles)) {
    indices <- .np_plreg_formula_indices(roles[[i]], joint)
    prediction <- predictions[[i]]$prediction
    for (j in seq_along(indices)) {
      index <- indices[[j]]
      if (seen[[index]] && !identical(combined[[index + 1L]], prediction[[j + 1L]]))
        stop("inconsistent partially linear trained prediction terms", call. = FALSE)
      combined[[index + 1L]] <- prediction[[j + 1L]]
      seen[[index]] <- TRUE
    }
  }
  if (!all(seen))
    stop("incomplete partially linear formula metadata", call. = FALSE)
  attr(joint, "predvars") <- combined
  joint
}

.np_plreg_formula_split <- function(frame, terms, xterms) {
  joint <- attr(frame, "terms")
  split.role <- function(role) {
    indices <- .np_plreg_formula_indices(role, joint)
    out <- frame[, indices, drop = FALSE]
    attr(role, "predvars") <- as.call(c(list(quote(list)),
      as.list(attr(joint, "predvars"))[indices + 1L]))
    classes <- attr(joint, "dataClasses")
    if (!is.null(classes)) attr(role, "dataClasses") <- classes[indices]
    attr(out, "terms") <- role
    attr(out, "na.action") <- attr(frame, "na.action")
    out
  }
  list(yz = split.role(terms), x = split.role(xterms))
}

.np_plreg_formula_frame <- function(bws, data = NULL) {
  m <- match(c("formula", "data", "subset", "na.action"),
             names(bws$call), nomatch = 0L)
  args <- as.list(bws$call[c(1L, m)])[-1L]
  args$formula <- .np_plreg_formula_terms(bws)
  if (!is.null(data)) args$data <- data
  .np_bws_formula_model_frame(bws, args, data.override = !is.null(data))
}

.np_plreg_formula_training <- function(bws, data = NULL) {
  frame <- .np_plreg_formula_frame(bws, data)
  roles <- .np_plreg_formula_split(frame, bws$terms, bws$xterms)
  list(txdat = roles$x, tydat = model.response(roles$yz),
       tzdat = roles$yz[, bws$chromoly[[3L]], drop = FALSE])
}
