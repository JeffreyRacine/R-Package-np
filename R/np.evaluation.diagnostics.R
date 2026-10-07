# Outcome discovery and diagnostic provenance for npindex and npconmode only.
# Evaluation outcomes never choose prediction rows or enter the fitted model.
.np_diagnostics_response <- function(newdata, tt = NULL, ynames = NULL,
                                     required = FALSE) {
  newdata <- toFrame(newdata)
  if (!is.null(tt)) {
    response <- attr(terms(tt), "variables")[[2L]]
    required.names <- all.vars(response)
  } else {
    required.names <- ynames
  }
  present <- length(required.names) > 0L &&
    all(required.names %in% names(newdata))
  if (!present && !required)
    return(NULL)
  if (!length(required.names))
    stop("evaluation response names are unavailable; supply 'eydat' explicitly",
         call. = FALSE)
  if (any(vapply(required.names, function(nm) sum(names(newdata) == nm) > 1L, logical(1))))
    stop("evaluation response columns must have unique names", call. = FALSE)
  npValidateNewdataColumns(newdata, required.names)
  if (!is.null(tt)) {
    # Validate above before evaluation can consult the formula environment.
    # Evaluate only the response: model.frame's raw-row warning is inappropriate
    # for length-changing time series which are aligned by time below.
    value <- eval(response, newdata, environment(tt))
  } else {
    if (length(required.names) != 1L)
      stop("evaluation response must have one unambiguous column", call. = FALSE)
    value <- newdata[[required.names]]
  }
  if (!inherits(value, "ts") && NROW(value) != NROW(newdata))
    stop("evaluation response must have one value per 'newdata' row", call. = FALSE)
  value
}

# Retain prediction time coordinates without evaluating RHS transformations twice.
# The shared frame owner still applies its existing alignment and NA policy.
.np_diagnostics_model_frame <- function(tt, newdata) {
  prediction <- .np_formula_unwrap_prediction(tt)$prediction
  values <- eval(prediction, newdata, environment(tt))
  series <- Filter(function(x) inherits(x, "ts"), values)
  times <- NULL
  if (length(series)) {
    grids <- lapply(series, function(x) {
      stats::ts(as.numeric(stats::time(x)), start = stats::start(x),
                frequency = stats::frequency(x))
    })
    grid <- if (length(grids) == 1L) grids[[1L]] else
      do.call(stats::ts.intersect, grids)
    times <- as.numeric(stats::time(grid))
  }
  attr(tt, "predvars") <- values
  frame <- .np_formula_model_frame(tt, data = newdata)
  retained <- attr(frame, "terms")
  attr(retained, "predvars") <- prediction
  attr(frame, "terms") <- retained
  attr(frame, ".np.diagnostics.times") <- times
  frame
}

.np_diagnostics_align_response <- function(value, frame) {
  times <- attr(frame, ".np.diagnostics.times", exact = TRUE)
  if (inherits(value, "ts") && length(times)) {
    # Scores have no observation beyond the supplied response's time support.
    position <- (times - as.numeric(stats::time(value))[1L]) * stats::frequency(value)
    index <- as.integer(round(position)) + 1L
    available <- abs(position - round(position)) < 1e-7 &
      index >= 1L & index <= NROW(value)
    aligned <- rep(NA_real_, length(times))
    aligned[available] <- as.numeric(value)[index[available]]
    value <- aligned
  }
  omitted <- attr(frame, "na.action")
  if (length(omitted)) value <- if (is.null(dim(value)))
    value[-as.integer(omitted)] else value[-as.integer(omitted), , drop = FALSE]
  if (NROW(value) != NROW(frame))
    stop("evaluation response and predictor rows are not aligned", call. = FALSE)
  value
}

.np_diagnostics_summary <- function(object) {
  sample <- object[["diagnostics.sample", exact = TRUE]]
  if (is.null(sample))
    return(invisible(NULL)) # Historical objects do not identify their score sample.
  if (identical(sample, "unavailable")) {
    cat("\nEvaluation diagnostics unavailable: no observed outcome/prediction pairs.\n")
  } else {
    cat("\n", if (identical(sample, "training")) "Training" else "Evaluation",
        " diagnostics (", object[["diagnostics.nobs", exact = TRUE]],
        " scored observations)\n", sep = "")
  }
  invisible(NULL)
}

# Only automatically discovered scoring data may fail without aborting a fit.
# Explicit outcomes retain their validation errors and warnings.
.np_diagnostics_optional <- function(expr, required = FALSE) {
  if (required) return(expr)
  tryCatch(suppressWarnings(expr), error = function(e) NULL)
}

.np_index_diagnostics_response <- function(value, method, required = FALSE) {
  if (is.null(value)) return(NULL)
  .np_diagnostics_optional({
    if (match.arg(method, c("ichimura", "kleinspady")) == "kleinspady")
      .npindex_check_binary_response(value, "npindex() evaluation response")
    else if (!(is.numeric(value) || is.factor(value)) || !is.null(dim(value)))
      stop("npindex() evaluation response must be a numeric vector or factor",
           call. = FALSE)
    value
  }, required = required)
}

.np_index_training_mse <- function(object) {
  value <- object[["training.MSE", exact = TRUE]]
  if (!is.null(value)) return(value)
  # Historical objects lack this provenance; preserve their stored scale.
  sample <- object[["diagnostics.sample", exact = TRUE]]
  if (is.null(sample) || identical(sample, "training")) object$MSE else NA_real_
}
