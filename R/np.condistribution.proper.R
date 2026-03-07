.np_condist_proper_shadow_enabled <- function() {
  isTRUE(getOption("np.condist.proper.shadow", FALSE))
}

.np_condist_proper_reason_message <- function(reason, where = "npcdist()") {
  switch(
    as.character(reason),
    not_requested = "proper-distribution repair was not requested",
    no_eval_grid = "proper=TRUE requires explicit evaluation data that define a full y-grid for each fixed x; paired row-wise evaluation is insufficient",
    y_not_univariate_continuous = "proper=TRUE currently supports only univariate continuous y",
    gradients_unsupported = "proper=TRUE is currently unsupported when gradients=TRUE",
    x_slices_not_repeated = "proper=TRUE requires repeated fixed-x slices in the evaluation data",
    y_grid_not_common = "proper=TRUE requires every fixed-x slice to share a common ordered y-grid",
    y_grid_not_ordered = "proper=TRUE requires at least two strictly increasing y-grid points in each fixed-x slice",
    invalid_quadrature = "proper=TRUE could not construct valid grid weights from the evaluation y-grid",
    repair_disabled_shadow_mode = sprintf("proper=TRUE repair is disabled in internal shadow-validation mode in %s", where),
    sprintf("proper=TRUE is unsupported in %s for reason '%s'", where, as.character(reason))
  )
}

.np_condist_normalize_proper_control <- function(proper.control) {
  if (is.null(proper.control))
    proper.control <- list()

  if (!is.list(proper.control))
    stop("'proper.control' must be a list")

  ctrl <- list(
    tol = 1e-10,
    grid.check = TRUE,
    store.raw = TRUE,
    fail.on.unsupported = FALSE
  )

  known <- names(ctrl)
  if (length(proper.control)) {
    extra <- setdiff(names(proper.control), known)
    if (length(extra))
      stop(sprintf("unsupported entries in 'proper.control': %s",
                   paste(shQuote(extra), collapse = ", ")))
    ctrl[names(proper.control)] <- proper.control
  }

  ctrl$tol <- npValidatePositiveFiniteNumeric(ctrl$tol, "proper.control$tol")
  ctrl$grid.check <- npValidateScalarLogical(ctrl$grid.check, "proper.control$grid.check")
  ctrl$store.raw <- npValidateScalarLogical(ctrl$store.raw, "proper.control$store.raw")
  ctrl$fail.on.unsupported <- npValidateScalarLogical(
    ctrl$fail.on.unsupported,
    "proper.control$fail.on.unsupported"
  )

  ctrl
}

.np_condist_validate_proper_args <- function(proper = FALSE,
                                             proper.method = c("isotonic"),
                                             proper.control = list()) {
  proper <- npValidateScalarLogical(proper, "proper")
  proper.method <- match.arg(as.character(proper.method)[1L], c("isotonic"))
  proper.control <- .np_condist_normalize_proper_control(proper.control)

  list(
    proper.requested = proper,
    proper.method = proper.method,
    proper.control = proper.control
  )
}

.np_condist_make_reason_info <- function(reason,
                                         supported = FALSE,
                                         slice.count = 0L,
                                         grid.common = FALSE,
                                         lower.bound = 0,
                                         upper.bound = 1,
                                         monotone.violations.raw = NULL,
                                         range.raw = NULL,
                                         projection.distance.l2 = NULL) {
  list(
    supported = supported,
    reason = as.character(reason)[1L],
    slice.count = as.integer(slice.count),
    grid.common = isTRUE(grid.common),
    lower.bound = lower.bound,
    upper.bound = upper.bound,
    monotone.violations.raw = monotone.violations.raw,
    range.raw = range.raw,
    projection.distance.l2 = projection.distance.l2
  )
}

.np_condist_weighted_pava <- function(f, w, tol = 1e-10) {
  f <- as.double(f)
  w <- as.double(w)

  if (length(f) != length(w))
    stop("length mismatch between distribution values and grid weights")
  if (any(!is.finite(f)))
    stop("distribution values must be finite")
  if (any(!is.finite(w)) || any(w <= 0))
    stop("grid weights must be finite and strictly positive")
  if (!is.finite(tol) || tol <= 0)
    stop("projection tolerance must be positive and finite")

  block.start <- seq_along(f)
  block.end <- seq_along(f)
  block.weight <- w
  block.value <- f
  nblock <- length(f)

  i <- 1L
  while (i < nblock) {
    if (block.value[i] <= block.value[i + 1L] + tol) {
      i <- i + 1L
      next
    }

    new.weight <- block.weight[i] + block.weight[i + 1L]
    new.value <- (block.weight[i] * block.value[i] +
                    block.weight[i + 1L] * block.value[i + 1L]) / new.weight

    block.end[i] <- block.end[i + 1L]
    block.weight[i] <- new.weight
    block.value[i] <- new.value

    if (i + 1L < nblock) {
      keep <- seq.int(i + 2L, nblock)
      block.start[(i + 1L):(nblock - 1L)] <- block.start[keep]
      block.end[(i + 1L):(nblock - 1L)] <- block.end[keep]
      block.weight[(i + 1L):(nblock - 1L)] <- block.weight[keep]
      block.value[(i + 1L):(nblock - 1L)] <- block.value[keep]
    }

    nblock <- nblock - 1L
    if (i > 1L)
      i <- i - 1L
  }

  fitted <- numeric(length(f))
  for (j in seq_len(nblock))
    fitted[block.start[j]:block.end[j]] <- block.value[j]

  fitted
}

.np_condist_project_bounded_isotonic <- function(f,
                                                 w,
                                                 lower = 0,
                                                 upper = 1,
                                                 tol = 1e-10) {
  f <- as.double(f)
  w <- as.double(w)

  if (!is.finite(lower) || !is.finite(upper) || lower > upper)
    stop("invalid lower/upper bounds for isotonic projection")

  big.weight <- max(1, sum(w)) * 1e12
  aug.fit <- .np_condist_weighted_pava(
    f = c(lower, f, upper),
    w = c(big.weight, w, big.weight),
    tol = tol
  )

  out <- aug.fit[-c(1L, length(aug.fit))]
  out <- pmin(pmax(out, lower), upper)

  if (any(diff(out) < -10 * tol))
    stop("bounded isotonic projection failed monotonicity check")

  out
}

.np_condist_apply_proper <- function(object,
                                     proper.method = "isotonic",
                                     proper.control = list()) {
  proper.method <- match.arg(as.character(proper.method)[1L], c("isotonic"))
  proper.control <- .np_condist_normalize_proper_control(proper.control)

  grid <- .np_condens_detect_proper_grid(object = object, tol = proper.control$tol)
  if (!isTRUE(grid$supported)) {
    info <- .np_condist_make_reason_info(
      reason = grid$reason,
      supported = FALSE,
      slice.count = 0L,
      grid.common = FALSE
    )
    return(list(
      applied = FALSE,
      reason = grid$reason,
      proper.info = info
    ))
  }

  raw <- as.double(object$condist)
  repaired <- raw
  monotone.violations.raw <- integer(length(grid$slices))
  range.raw <- matrix(NA_real_, nrow = length(grid$slices), ncol = 2L,
                      dimnames = list(NULL, c("min", "max")))
  projection.distance <- numeric(length(grid$slices))

  for (i in seq_along(grid$slices)) {
    idx <- grid$slices[[i]]
    f.slice <- raw[idx]
    w.slice <- tryCatch(
      .np_condens_trapezoid_weights(grid$y.grid, tol = proper.control$tol),
      error = function(e) e
    )
    if (inherits(w.slice, "error")) {
      info <- .np_condist_make_reason_info(
        reason = "invalid_quadrature",
        supported = FALSE,
        slice.count = length(grid$slices),
        grid.common = TRUE
      )
      return(list(
        applied = FALSE,
        reason = "invalid_quadrature",
        proper.info = info
      ))
    }

    g.slice <- .np_condist_project_bounded_isotonic(
      f = f.slice,
      w = w.slice,
      lower = 0,
      upper = 1,
      tol = proper.control$tol
    )

    repaired[idx] <- g.slice
    monotone.violations.raw[i] <- sum(diff(f.slice) < -proper.control$tol)
    range.raw[i, ] <- c(min(f.slice), max(f.slice))
    projection.distance[i] <- sqrt(sum(w.slice * (g.slice - f.slice)^2))
  }

  info <- .np_condist_make_reason_info(
    reason = "applied",
    supported = TRUE,
    slice.count = length(grid$slices),
    grid.common = TRUE,
    monotone.violations.raw = monotone.violations.raw,
    range.raw = range.raw,
    projection.distance.l2 = projection.distance
  )

  list(
    applied = TRUE,
    condist = repaired,
    condist.raw = if (isTRUE(proper.control$store.raw)) raw else NULL,
    proper.info = info
  )
}

.np_condist_finalize_proper_object <- function(object,
                                               proper = FALSE,
                                               proper.method = c("isotonic"),
                                               proper.control = list(),
                                               where = "npcdist()") {
  args <- .np_condist_validate_proper_args(
    proper = proper,
    proper.method = proper.method,
    proper.control = proper.control
  )

  object$proper.requested <- args$proper.requested
  object$proper.applied <- FALSE
  object$proper.method <- if (isTRUE(args$proper.requested)) args$proper.method else NULL
  object$condist.raw <- NULL
  object$proper.info <- .np_condist_make_reason_info(
    reason = if (isTRUE(args$proper.requested)) "pending" else "not_requested",
    supported = FALSE
  )

  if (!isTRUE(args$proper.requested))
    return(object)

  if (.np_condist_proper_shadow_enabled()) {
    object$proper.info <- .np_condist_make_reason_info(
      reason = "repair_disabled_shadow_mode",
      supported = FALSE
    )
    return(object)
  }

  proper.out <- .np_condist_apply_proper(
    object = object,
    proper.method = args$proper.method,
    proper.control = args$proper.control
  )

  if (!isTRUE(proper.out$applied)) {
    object$proper.info <- proper.out$proper.info
    if (isTRUE(args$proper.control$fail.on.unsupported)) {
      stop(.np_condist_proper_reason_message(
        reason = proper.out$reason,
        where = where
      ), call. = FALSE)
    }
    return(object)
  }

  object$condist <- proper.out$condist
  object$proper.applied <- TRUE
  object$proper.info <- proper.out$proper.info
  if (!is.null(proper.out$condist.raw))
    object$condist.raw <- proper.out$condist.raw

  object
}
