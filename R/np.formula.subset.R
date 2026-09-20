# terms() quotes non-syntactic symbols; model.frame() uses their actual names.
# Resolve only a simple symbol. Never evaluate a term or rewrite a transform,
# interaction, formula, predvars expression, or the user's data-frame names.
.np_formula_term_names <- function(labels) {
  quoted <- which(startsWith(labels, "`"))
  if (!length(quoted)) return(labels)
  for (i in quoted) {
    term <- str2lang(labels[[i]])
    if (is.symbol(term)) labels[[i]] <- as.character(term)
  }
  labels
}

# Inspect terms, not text: offset as a variable/response or inside an ordinary
# transformation is not an offset special and retains its existing semantics.
.np_formula_validate_terms <- function(tt) {
  if (length(attr(tt, "offset")))
    stop("offset() terms are not supported in this estimator formula; remove the offset term or use explicitly prepared data",
         call. = FALSE)
  tt
}

# Call only after the existing owner has resolved its formula promise. Never
# evaluate a construction call or a response/covariate expression here.
.np_formula_validate_syntax <- function(formula, conditional.response = FALSE) {
  if (conditional.response && length(formula) == 3L) {
    named.response <- function(expr) {
      if (is.symbol(expr)) return(TRUE)
      if (!is.call(expr)) return(FALSE)
      if (identical(expr[[1L]], quote(`(`)) && length(expr) == 2L)
        return(named.response(expr[[2L]]))
      if (identical(expr[[1L]], quote(`+`)) && length(expr) %in% c(2L, 3L))
        return(all(vapply(as.list(expr)[-1L], named.response, logical(1L))))
      FALSE
    }
    if (!named.response(formula[[2L]]))
      stop("conditional formula responses must be variable names separated by '+'; create an explicit transformed variable in data or use the native data interface",
           call. = FALSE)
  }
  parts <- function(expr) {
    if (is.call(expr) && identical(expr[[1L]], quote(`|`)) &&
        length(expr) == 3L)
      return(c(parts(expr[[2L]]), parts(expr[[3L]])))
    list(expr)
  }
  for (rhs in parts(formula[[length(formula)]])) {
    role <- as.call(list(quote(`~`), rhs))
    class(role) <- "formula"
    environment(role) <- environment(formula)
    .np_formula_validate_terms(terms(role, allowDotAsName = TRUE))
  }
  invisible(formula)
}

# Formula-owned subset is non-standard evaluation: dispatch must inspect other
# arguments without forcing it outside the data mask. Read retained arguments
# from their original promises, never by re-evaluating substituted expressions.
.np_formula_subset_indices <- function(method, expressions) {
  dot.names <- names(expressions)
  indices <- if (is.null(dot.names)) integer() else
    which(!is.na(pmatch(dot.names, "subset")))
  if (!is.null(method) && length(expressions)) {
    markers <- lapply(seq_along(expressions), function(i)
      as.name(paste0(".np_formula_dispatch_dot_", i)))
    names(markers) <- dot.names
    synthetic <- as.call(c(list(as.name(".np_formula_dispatch")),
      list(bws = as.name(".np_formula_dispatch_bws")), markers))
    matched <- tryCatch(match.call(definition = method, call = synthetic,
                                  expand.dots = FALSE),
                        error = function(e) NULL)
    if (!is.null(matched) && "subset" %in% names(matched))
      indices <- union(indices, which(vapply(markers, identical, logical(1L),
                                             matched[["subset"]])))
  }
  indices
}

.np_formula_dispatch_args <- function(method, expressions, promise.frame) {
  omitted <- .np_formula_subset_indices(method, expressions)
  if (!length(omitted))
    return(eval(quote(list(...)), envir = promise.frame))
  indices <- setdiff(seq_along(expressions), omitted)
  args <- lapply(indices, function(i)
    eval(substitute(...elt(index), list(index = i)), envir = promise.frame))
  if (!is.null(names(expressions)))
    names(args) <- names(expressions)[indices]
  args
}

# Existing value-list constructor calls still need the original subset syntax.
# The constructor's model.frame, not this helper, evaluates that expression.
.np_formula_dispatch_call <- function(fun, args, expressions, envir) {
  indices <- .np_formula_subset_indices(NULL, expressions)
  do.call(fun, c(args, as.list(expressions)[indices]), envir = envir)
}

# Evaluate expression values once, align them transiently, and retain portable
# prediction expressions (including makepredictcall metadata), never data in
# those prediction expressions. Bandwidth owners separately retain the frame.
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
                                    drop.unused.levels = FALSE, xlev = NULL, ...,
                                    .np.capture = NULL, .np.auxiliary = list()) {
  if (!is.data.frame(data) && !is.environment(data) && !is.null(attr(data, "class")))
    data <- as.data.frame(data)
  tt <- if (inherits(formula, "terms")) formula else terms(formula, data = data)
  .np_formula_validate_terms(tt)
  prepared <- .np_formula_unwrap_prediction(tt)
  prediction <- prepared$prediction
  values <- eval(prediction, data, environment(tt))
  if (prepared$makepredictcall) {
    variables <- attr(tt, "variables")
    for (i in seq_along(values))
      prediction[[i + 1L]] <- stats::makepredictcall(values[[i]], variables[[i + 1L]])
  }
  # Observation-level auxiliaries share the formula's time intersection,
  # subset and NA action. They are model-frame extras, never predictors.
  aligned <- .np_formula_align_values(c(values, .np.auxiliary))
  # I(numeric_vector) protects formula syntax, not a distinct kernel data type.
  # Keep matrices and all other explicit classes subject to existing policy.
  for (i in seq_along(values)) {
    value <- aligned[[i]]
    if (identical(attr(value, "class"), "AsIs") &&
        is.numeric(value) && is.null(dim(value))) {
      class(value) <- NULL
      aligned[[i]] <- value
    }
  }
  attr(tt, "predvars") <- aligned[seq_along(values)]
  auxiliary <- aligned[length(values) + seq_along(.np.auxiliary)]
  selected <- if (length(auxiliary))
    vapply(auxiliary, NROW, integer(1L)) != NROW(aligned[[1L]]) else logical()
  frame.call <- match.call()
  frame.call[[1L]] <- quote(stats::model.frame)
  frame.call[["formula"]] <- tt
  frame.call["data"] <- list(data)
  frame.call$.np.capture <- NULL
  frame.call$.np.auxiliary <- NULL
  if (length(auxiliary))
    frame.call[names(auxiliary)[!selected]] <- auxiliary[!selected]
  if (!is.null(.np.capture)) {
    # Freeze the resolved policy, not its caller-side expression. Force the
    # promise once, and let model.frame invoke the policy exactly once.
    policy <- if (!missing(na.action)) na.action else {
      action <- attr(data, "na.action")
      if (!is.null(action) && mode(action) != "numeric") action else
        getOption("na.action", stats::na.fail)
    }
    if (is.character(policy))
      policy <- get(policy, envir = asNamespace("stats"), mode = "function")
    .np.capture$na.action <- policy
    frame.call["na.action"] <- list(policy)
  }
  frame <- eval(frame.call, parent.frame())
  # Preserve unambiguously preselected legacy auxiliaries. Full-length inputs
  # always belong to the original sample, even when subset only reorders it.
  for (name in names(auxiliary)[selected]) {
    value <- auxiliary[[name]]
    if (NROW(value) != nrow(frame))
      stop(sprintf("'%s' must match the original or selected formula sample", name),
           call. = FALSE)
    frame[[paste0("(", name, ")")]] <- value
  }
  retained <- attr(frame, "terms")
  attr(retained, "predvars") <- prediction
  attr(frame, "terms") <- retained
  if (!is.null(.np.capture)) frame <- .np_formula_complete_training_frame(frame)
  frame
}

# Native training owners use complete cases even when a model-frame policy
# leaves missing values. Record that existing exclusion before formula owners
# freeze the sample, rather than overwriting native omission bookkeeping later.
.np_formula_complete_training_frame <- function(frame) {
  missing.rows <- which(!stats::complete.cases(frame))
  if (!length(missing.rows)) return(frame)
  previous <- attr(frame, "na.action")
  kept <- seq_len(nrow(frame) + length(previous))
  if (length(previous)) kept <- kept[-as.integer(previous)]
  omitted <- c(as.integer(previous), kept[missing.rows])
  labels <- c(names(previous), row.names(frame)[missing.rows])
  order <- order(omitted)
  action <- structure(omitted[order], names = labels[order],
    class = if (inherits(previous, "exclude")) "exclude" else "omit")
  result <- frame[-missing.rows, , drop = FALSE]
  attr(result, "terms") <- attr(frame, "terms")
  attr(result, "na.action") <- action
  result
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

.np_bws_retain_formula_training <- function(bws, frame, na.action) {
  attr(frame, ".np.na.policy") <- NULL
  bws[[".np.formula.training"]] <- list(frame = frame, na.action = na.action)
  bws[[".np.native.training"]] <- NULL
  call.env <- environment(bws$call)
  # An internal constructor activation is not the owner of a user's formula.
  # Retained training values make this transient frame unnecessary; preserve
  # the original formula environment for transforms and explicit-data refits.
  environment(bws$call) <- .np_call_owner_environment(
    call.env, replacement = environment(bws$formula))
  bws
}

# Project stored formula metadata, never formula expressions, into native roles.
.np_bws_formula_roles <- function(bws) {
  tt <- bws[["terms", exact = TRUE]]
  if (inherits(bws, c("bandwidth", "dbandwidth")))
    return(list(dat = .np_formula_term_names(attr(tt, "term.labels"))))
  if (inherits(bws, c("conbandwidth", "condbandwidth")))
    return(list(xdat = bws$variableNames[["terms"]],
                ydat = bws$variableNames[["response"]]))
  # A transformed scalar response has its full term label, not all.vars().
  variable <- attr(tt, "variables")[[2L]]
  response <- if (is.symbol(variable)) as.character(variable) else
    paste(deparse(variable, width.cutoff = 500L), collapse = "")
  chromoly <- bws[["chromoly", exact = TRUE]]
  if (inherits(bws, c("plbandwidth", "scbandwidth"))) {
    out <- list(xdat = .np_formula_term_names(chromoly[[2L]]), ydat = response)
    if (length(chromoly) == 3L) out$zdat <- .np_formula_term_names(chromoly[[3L]])
    return(out)
  }
  list(xdat = .np_formula_term_names(attr(tt, "term.labels")), ydat = response)
}

.np_bws_retain_formula_roles <- function(bws, training) {
  roles <- .np_bws_formula_roles(bws)
  if (!setequal(names(roles), names(training)))
    stop("incomplete retained formula training roles", call. = FALSE)
  retained <- bws[[".np.formula.training", exact = TRUE]]
  old <- retained[["frame", exact = TRUE]]
  columns <- list()
  rows <- NULL
  for (role in names(roles)) {
    values <- toFrame(training[[role]])
    if (ncol(values) != length(roles[[role]]))
      stop("retained formula training role has incompatible columns", call. = FALSE)
    if (is.null(rows)) rows <- row.names(values)
    if (nrow(values) != length(rows))
      stop("retained formula training roles have incompatible rows", call. = FALSE)
    for (j in seq_along(roles[[role]])) {
      name <- roles[[role]][[j]]
      if (name %in% names(columns) && !identical(columns[[name]], values[[j]]))
        stop("overlapping formula training roles disagree", call. = FALSE)
      columns[[name]] <- values[[j]]
    }
  }
  if (!is.null(old) && setequal(names(old), names(columns)) &&
      identical(row.names(old), rows) &&
      all(vapply(names(old), function(name) identical(old[[name]], columns[[name]]), logical(1L)))) {
    bws[[".np.native.training"]] <- NULL
    return(bws)
  }
  tt <- if (inherits(bws, "plbandwidth")) .np_plreg_formula_terms(bws) else bws$terms
  column.order <- if (!is.null(old)) names(old) else
    vapply(as.list(attr(tt, "variables"))[-1L], function(x)
      paste(deparse(x, width.cutoff = 500L), collapse = ""), character(1L))
  if (!setequal(column.order, names(columns)))
    stop("retained formula training columns do not match terms", call. = FALSE)
  columns <- columns[column.order]
  frame <- as.data.frame(columns, optional = TRUE, row.names = rows)
  # Native replacement inputs already represent the formula's evaluated terms.
  # Retain their complete cases without replaying old subset/NA expressions.
  attr(frame, "terms") <- tt
  frame <- stats::na.omit(frame)
  policy <- if (is.null(retained)) stats::na.omit else retained[["na.action", exact = TRUE]]
  .np_bws_retain_formula_training(bws, frame, policy)
}

# A partial response replacement conditions on the already-selected design.
# Its rows therefore index the retained frame, never the original data before
# subset/NA selection. Do not replay expressions or infer ownership from length.
.np_formula_replace_response <- function(frame, bws, response, overrides,
                                         data.override = FALSE,
                                         matrix.response = FALSE) {
  if (data.override)
    stop("partial response replacement cannot be combined with 'data'; put the replacement response in 'data' instead",
         call. = FALSE)
  roles <- .np_bws_formula_roles(bws)
  columns <- roles[["ydat"]]
  if (!length(columns) || !all(columns %in% names(frame)))
    stop("partial response replacement requires retained response columns", call. = FALSE)
  if (length(intersect(columns, unlist(roles[names(roles) != "ydat"], use.names = FALSE))))
    stop("a response used as a predictor requires full 'data' replacement", call. = FALSE)
  if (NROW(response) != nrow(frame))
    stop("replacement response must have one row per retained training observation, in retained-sample order; use 'data' to replace the original sample",
         call. = FALSE)
  if (isTRUE(matrix.response) && length(columns) == 1L && is.matrix(response)) {
    frame[[columns]] <- I(response)
  } else {
    values <- toFrame(response)
    if (ncol(values) != length(columns))
      stop("replacement response has an incompatible number of columns", call. = FALSE)
    for (j in seq_along(columns)) frame[[columns[[j]]]] <- values[[j]]
  }
  retained <- bws[[".np.formula.training", exact = TRUE]]
  policy <- if ("na.action" %in% names(overrides)) overrides[["na.action"]] else
    if (!is.null(retained)) retained[["na.action", exact = TRUE]] else {
      effective <- attr(frame, ".np.na.policy", exact = TRUE)
      if (!is.null(effective)) effective[[1L]] else getOption("na.action", stats::na.omit)
    }
  if (is.character(policy))
    policy <- get(policy, envir = asNamespace("stats"), mode = "function")
  tt <- attr(frame, "terms")
  attr(frame, "na.action") <- NULL
  if (!is.null(policy)) frame <- policy(frame)
  attr(frame, "terms") <- tt
  frame <- .np_formula_complete_training_frame(frame)
  attr(frame, ".np.na.policy") <- list(policy)
  frame
}

.np_bws_retain_fit_frame <- function(bws, frame) {
  retained <- bws[[".np.formula.training", exact = TRUE]]
  effective <- attr(frame, ".np.na.policy", exact = TRUE)
  policy <- if (!is.null(effective)) effective[[1L]] else if (is.null(retained)) getOption("na.action", stats::na.omit) else
    retained[["na.action", exact = TRUE]]
  .np_bws_retain_formula_training(bws, frame, policy)
}

.np_formula_default_call <- function(call, definition, caller, required.training = NULL) {
  # Name the bandwidth before changing the callee. Keep the remaining call
  # shape: an unnamed formula beside bws= belongs to the formula dispatcher,
  # not to an explicitly named native training argument.
  call <- .np_formula_expand_call(call, caller)
  matched <- match.call(definition = definition, call = call, expand.dots = TRUE)
  if (!"bws" %in% names(matched) || "bws" %in% names(call)) return(call)
  labels <- names(call)
  if (is.null(labels)) labels <- rep.int("", length(call))
  labels[is.na(labels)] <- ""
  partial <- which(nzchar(labels) & startsWith("bws", labels))
  if (!length(partial)) {
    # These defaults also accept data-first calls (e.g. npudens(x) or
    # npreg(x, y)). Matching the fitter's first formal does not by itself
    # establish bandwidth ownership: its selector has a data-first NULL
    # method. Preserve that syntax unless the training slots establish bws.
    training <- if (is.null(required.training))
      grep("^t(x|y|z)?dat$", names(formals(definition)), value = TRUE) else required.training
    first <- training[[1L]]
    first.named <- any(nzchar(labels) &
      (startsWith(first, labels) | labels == substring(first, 2L)))
    complete <- all(training %in% names(matched))
    formula.after <- !first.named && !complete && first %in% names(matched) &&
      inherits(get(first, envir = parent.frame(), inherits = FALSE), "formula")
    if (!first.named && !complete && !formula.after && !"formula" %in% labels)
      return(call)
  }
  index <- if (length(partial)) partial[[1L]] else which(labels[-1L] == "")[[1L]] + 1L
  labels[[index]] <- "bws"
  names(call) <- labels
  call
}

.np_bws_formula_model_frame <- function(bws, mf.args, data.override = FALSE,
                                        overrides = list()) {
  if (inherits(bws, c("conbandwidth", "condbandwidth")))
    .np_formula_validate_syntax(bws$formula, conditional.response = TRUE)
  training <- bws[[".np.formula.training", exact = TRUE]]
  if (!is.null(training)) {
    if (!is.list(training) || !is.data.frame(training[["frame", exact = TRUE]]) ||
        !("na.action" %in% names(training)))
      stop("invalid retained formula training state", call. = FALSE)
    if (!data.override)
      return(.np_formula_complete_training_frame(training[["frame", exact = TRUE]]))
    mf.args["na.action"] <- training["na.action"]
  }
  # Compatibility for objects saved before training frames were retained.
  # Their original values cannot be reconstructed after caller rebinding.
  call.env <- environment(bws$call)
  if (is.null(training) && is.environment(call.env)) {
    # These are the value arguments extracted from the saved model-frame call.
    # Resolve them lazily in their original owner, including forwarded ..n
    # promises. Formula variables and subset still use the formula/data mask.
    # eval is necessary because do.call below runs in the formula environment,
    # which need not contain the constructor's data or NA-action bindings.
    owned <- intersect(c(if (!data.override) "data", "na.action"), names(mf.args))
    for (name in owned) {
      if (is.language(mf.args[[name]]))
        mf.args[[name]] <- substitute(base::eval(quote(EXPR), envir = OWNER),
                                     list(EXPR = mf.args[[name]], OWNER = call.env))
    }
  }
  if (data.override && "na.action" %in% names(overrides))
    mf.args["na.action"] <- overrides["na.action"]
  capture <- new.env(parent = emptyenv())
  mf.args$.np.capture <- capture
  frame <- do.call(.np_formula_model_frame, mf.args,
                  envir = environment(mf.args[["formula"]]))
  # Invocation-local metadata carries even an explicit NULL policy. The
  # retention owner strips this attribute from the stored model frame.
  attr(frame, ".np.na.policy") <- list(capture$na.action)
  frame
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

.np_plreg_formula_frame <- function(bws, data = NULL, overrides = list()) {
  m <- match(c("formula", "data", "subset", "na.action"),
             names(bws$call), nomatch = 0L)
  args <- as.list(bws$call[c(1L, m)])[-1L]
  args$formula <- .np_plreg_formula_terms(bws)
  if (!is.null(data)) args$data <- data
  .np_bws_formula_model_frame(bws, args, data.override = !is.null(data), overrides = overrides)
}

.np_plreg_formula_training <- function(bws, data = NULL) {
  frame <- .np_plreg_formula_frame(bws, data)
  roles <- .np_plreg_formula_split(frame, bws$terms, bws$xterms)
  list(txdat = roles$x, tydat = model.response(roles$yz),
       tzdat = roles$yz[, .np_formula_term_names(bws$chromoly[[3L]]), drop = FALSE])
}
