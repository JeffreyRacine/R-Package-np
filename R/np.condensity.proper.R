.np_condens_proper_shadow_enabled <- function() {
  isTRUE(getOption("np.condens.proper.shadow", FALSE))
}

.np_condens_proper_reason_message <- function(reason, where = "npcdens()") {
  switch(
    as.character(reason),
    not_requested = "proper-density repair was not requested",
    already_proper = "proper=TRUE was requested, but the estimator is already proper by construction",
    no_eval_grid = "proper=TRUE requires explicit evaluation data that define a full y-grid for each fixed x; paired row-wise evaluation is insufficient",
    scope_not_selected = "proper=TRUE was requested, but 'proper.control$apply' does not target the values returned by this object",
    slice_disabled = "proper=TRUE slice repair is disabled by internal dispatcher controls",
    slice_context_missing = "proper=TRUE slice repair requires explicit evaluation and training data context",
    slice_invalid_master_grid = "proper=TRUE slice repair could not construct a valid internal y-grid",
    slice_eval_failed = "proper=TRUE slice repair failed while evaluating the internal explicit-grid oracle",
    y_not_univariate_continuous = "proper=TRUE currently supports only univariate continuous y",
    gradients_unsupported = "proper=TRUE is currently unsupported when gradients=TRUE",
    x_slices_not_repeated = "proper=TRUE requires repeated fixed-x slices in the evaluation data",
    y_grid_not_common = "proper=TRUE requires every fixed-x slice to share a common ordered y-grid",
    y_grid_not_ordered = "proper=TRUE requires at least two strictly increasing y-grid points in each fixed-x slice",
    invalid_quadrature = "proper=TRUE could not construct valid trapezoidal quadrature weights from the evaluation y-grid",
    repair_disabled_shadow_mode = sprintf("proper=TRUE repair is disabled in internal shadow-validation mode in %s", where),
    sprintf("proper=TRUE is unsupported in %s for reason '%s'", where, as.character(reason))
  )
}

.np_condens_normalize_proper_control <- function(proper.control) {
  if (is.null(proper.control))
    proper.control <- list()

  if (!is.list(proper.control))
    stop("'proper.control' must be a list")

  ctrl <- list(
    tol = 1e-10,
    grid.check = TRUE,
    store.raw = TRUE,
    fail.on.unsupported = FALSE,
    mode = "grid",
    apply = "evaluation",
    slice.grid.size = 101L,
    slice.extend.factor = 0.1
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
  ctrl$mode <- match.arg(as.character(ctrl$mode)[1L], c("grid", "slice"))
  ctrl$apply <- match.arg(as.character(ctrl$apply)[1L], c("evaluation", "fitted", "both"))
  ctrl$`slice.grid.size` <- npValidatePositiveInteger(
    ctrl$`slice.grid.size`,
    "proper.control$slice.grid.size"
  )
  ctrl$`slice.extend.factor` <- .np_condens_validate_nonnegative_finite_numeric(
    ctrl$`slice.extend.factor`,
    "proper.control$slice.extend.factor"
  )

  ctrl
}

.np_condens_validate_proper_args <- function(proper = FALSE,
                                             proper.method = c("project"),
                                             proper.control = list()) {
  proper <- npValidateScalarLogical(proper, "proper")
  proper.method <- match.arg(as.character(proper.method)[1L], c("project"))
  proper.control <- .np_condens_normalize_proper_control(proper.control)

  list(
    proper.requested = proper,
    proper.method = proper.method,
    proper.control = proper.control
  )
}

.np_condens_make_reason_info <- function(reason,
                                         supported = FALSE,
                                         slice.count = 0L,
                                         grid.common = FALSE,
                                         quadrature = "trapezoid",
                                         negative.count.raw = NULL,
                                         integral.raw = NULL,
                                         projection.distance.l2 = NULL) {
  list(
    supported = supported,
    reason = as.character(reason)[1L],
    slice.count = as.integer(slice.count),
    grid.common = isTRUE(grid.common),
    quadrature = quadrature,
    negative.count.raw = negative.count.raw,
    integral.raw = integral.raw,
    projection.distance.l2 = projection.distance.l2
  )
}

.np_condens_trapezoid_weights <- function(y, tol = 1e-10) {
  y <- as.double(y)
  if (length(y) < 2L)
    stop("need at least two y-grid points to construct trapezoidal weights")
  if (any(!is.finite(y)))
    stop("y-grid points must be finite")

  dy <- diff(y)
  if (any(!is.finite(dy)) || any(dy <= tol))
    stop("y-grid points must be strictly increasing")

  w <- numeric(length(y))
  w[1L] <- dy[1L] / 2
  w[length(y)] <- dy[length(dy)] / 2
  if (length(y) > 2L)
    w[2:(length(y) - 1L)] <- (dy[-length(dy)] + dy[-1L]) / 2

  if (any(!is.finite(w)) || any(w <= 0))
    stop("invalid trapezoidal quadrature weights")

  w
}

.np_condens_project_weighted_simplex <- function(f,
                                                 w,
                                                 mass = 1,
                                                 tol = 1e-10,
                                                 maxit = 200L) {
  f <- as.double(f)
  w <- as.double(w)
  mass <- as.double(mass)[1L]
  maxit <- as.integer(maxit)[1L]

  if (length(f) != length(w))
    stop("length mismatch between density values and quadrature weights")
  if (any(!is.finite(f)))
    stop("density values must be finite")
  if (any(!is.finite(w)) || any(w <= 0))
    stop("quadrature weights must be finite and strictly positive")
  if (!is.finite(mass) || mass <= 0)
    stop("projection mass must be positive and finite")
  if (!is.finite(tol) || tol <= 0)
    stop("projection tolerance must be positive and finite")
  if (!is.finite(maxit) || maxit < 1L)
    stop("projection iteration count must be a positive integer")

  active.mass <- function(lambda) sum(w * pmax(f - lambda, 0))

  upper <- max(f)
  lower <- min(f) - max(1, max(abs(f)), mass / sum(w))

  mass.lower <- active.mass(lower)
  iter <- 0L
  while (mass.lower < mass && iter < maxit) {
    span <- max(1, abs(lower))
    lower <- lower - span
    mass.lower <- active.mass(lower)
    iter <- iter + 1L
  }

  if (mass.lower < mass)
    stop("could not bracket weighted simplex projection root")

  for (i in seq_len(maxit)) {
    mid <- (lower + upper) / 2
    mass.mid <- active.mass(mid)
    if (abs(mass.mid - mass) <= tol)
      break
    if (mass.mid > mass) {
      lower <- mid
    } else {
      upper <- mid
    }
  }

  lambda <- (lower + upper) / 2
  g <- pmax(f - lambda, 0)
  mass.g <- sum(w * g)

  if (!is.finite(mass.g) || abs(mass.g - mass) > 10 * tol)
    stop("weighted simplex projection failed to reach target mass")

  g
}

.np_condens_detect_proper_grid <- function(object, tol = 1e-10) {
  if (isTRUE(object$gradients)) {
    return(list(
      supported = FALSE,
      reason = "gradients_unsupported"
    ))
  }

  if (isTRUE(object$trainiseval)) {
    return(list(
      supported = FALSE,
      reason = "no_eval_grid"
    ))
  }

  if (!(identical(object$yndim, 1L) &&
        identical(object$yncon, 1L) &&
        identical(object$ynord, 0L) &&
        identical(object$ynuno, 0L))) {
    return(list(
      supported = FALSE,
      reason = "y_not_univariate_continuous"
    ))
  }

  xeval <- toFrame(object$xeval)
  yeval <- toFrame(object$yeval)

  if (nrow(yeval) < 2L) {
    return(list(
      supported = FALSE,
      reason = "y_grid_not_ordered"
    ))
  }

  if (ncol(xeval) == 0L) {
    groups <- factor(rep.int("all", nrow(yeval)))
  } else {
    groups <- do.call(interaction, c(unname(xeval), list(drop = TRUE, lex.order = TRUE)))
  }

  slices.raw <- split(seq_len(nrow(yeval)), groups, drop = TRUE)
  if (!length(slices.raw) || all(lengths(slices.raw) < 2L)) {
    return(list(
      supported = FALSE,
      reason = "x_slices_not_repeated"
    ))
  }

  yref <- NULL
  slices <- vector("list", length(slices.raw))

  for (i in seq_along(slices.raw)) {
    idx <- slices.raw[[i]]
    if (length(idx) < 2L) {
      return(list(
        supported = FALSE,
        reason = "x_slices_not_repeated"
      ))
    }

    y.slice <- as.double(yeval[[1L]][idx])
    ord <- order(y.slice)
    y.sorted <- y.slice[ord]

    if (length(y.sorted) < 2L || any(!is.finite(y.sorted)) ||
        any(diff(y.sorted) <= tol)) {
      return(list(
        supported = FALSE,
        reason = "y_grid_not_ordered"
      ))
    }

    if (is.null(yref)) {
      yref <- y.sorted
    } else if (!identical(y.sorted, yref)) {
      return(list(
        supported = FALSE,
        reason = "y_grid_not_common"
      ))
    }

    slices[[i]] <- idx[ord]
  }

  list(
    supported = TRUE,
    reason = NA_character_,
    slices = slices,
    y.grid = yref,
    slice.count = length(slices),
    grid.common = TRUE
  )
}

.np_condens_prepare_proper_plan <- function(object, proper.control = list()) {
  proper.control <- .np_condens_normalize_proper_control(proper.control)

  grid <- .np_condens_detect_proper_grid(object = object, tol = proper.control$tol)
  if (!isTRUE(grid$supported)) {
    info <- .np_condens_make_reason_info(
      reason = grid$reason,
      supported = FALSE,
      slice.count = 0L,
      grid.common = FALSE
    )
    return(list(
      supported = FALSE,
      reason = grid$reason,
      proper.info = info
    ))
  }

  weights <- tryCatch(
    .np_condens_trapezoid_weights(grid$y.grid, tol = proper.control$tol),
    error = function(e) e
  )
  if (inherits(weights, "error")) {
    info <- .np_condens_make_reason_info(
      reason = "invalid_quadrature",
      supported = FALSE,
      slice.count = length(grid$slices),
      grid.common = TRUE
    )
    return(list(
      supported = FALSE,
      reason = "invalid_quadrature",
      proper.info = info
    ))
  }

  list(
    supported = TRUE,
    slices = grid$slices,
    weights = weights,
    y.grid = grid$y.grid,
    proper.control = proper.control
  )
}

.np_condens_project_values_with_plan <- function(values,
                                                 plan,
                                                 progress.label = NULL) {
  if (!isTRUE(plan$supported))
    stop("proper projection plan is not supported")

  is.vector.input <- is.null(dim(values))
  values.mat <- if (is.vector.input) {
    matrix(as.double(values), nrow = 1L)
  } else {
    data.matrix(values)
  }

  if (ncol(values.mat) != sum(lengths(plan$slices)))
    stop("value length mismatch for proper density projection")

  out <- values.mat
  progress <- NULL
  if (!is.vector.input && !is.null(progress.label) && nrow(values.mat) > 1L) {
    progress <- .np_plot_stage_progress_begin(
      total = nrow(values.mat),
      label = as.character(progress.label)[1L]
    )
    on.exit(.np_plot_progress_end(progress), add = TRUE)
  }
  for (row in seq_len(nrow(values.mat))) {
    for (idx in plan$slices) {
      out[row, idx] <- .np_condens_project_weighted_simplex(
        f = values.mat[row, idx],
        w = plan$weights,
        mass = 1,
        tol = plan$proper.control$tol
      )
    }
    if (!is.null(progress))
      progress <- .np_plot_progress_tick(state = progress, done = row)
  }

  if (is.vector.input) as.vector(out[1L, ]) else out
}

.np_condens_apply_proper_grid <- function(object,
                                          proper.method = "project",
                                          proper.control = list()) {
  proper.method <- match.arg(as.character(proper.method)[1L], c("project"))
  proper.control <- .np_condens_normalize_proper_control(proper.control)

  plan <- .np_condens_prepare_proper_plan(
    object = object,
    proper.control = proper.control
  )
  if (!isTRUE(plan$supported)) {
    return(list(
      applied = FALSE,
      reason = plan$reason,
      proper.info = plan$proper.info
    ))
  }

  raw <- as.double(object$condens)
  repaired <- .np_condens_project_values_with_plan(raw, plan)
  negative.count.raw <- integer(length(plan$slices))
  integral.raw <- numeric(length(plan$slices))
  projection.distance <- numeric(length(plan$slices))

  for (i in seq_along(plan$slices)) {
    idx <- plan$slices[[i]]
    f.slice <- raw[idx]
    g.slice <- repaired[idx]
    negative.count.raw[i] <- sum(f.slice < 0)
    integral.raw[i] <- sum(plan$weights * f.slice)
    projection.distance[i] <- sqrt(sum(plan$weights * (g.slice - f.slice)^2))
  }

  info <- .np_condens_make_reason_info(
    reason = "applied",
    supported = TRUE,
    slice.count = length(plan$slices),
    grid.common = TRUE,
    negative.count.raw = negative.count.raw,
    integral.raw = integral.raw,
    projection.distance.l2 = projection.distance
  )

  list(
    applied = TRUE,
    condens = repaired,
    condens.raw = if (isTRUE(proper.control$store.raw)) raw else NULL,
    proper.info = info
  )
}

.np_condens_apply_proper <- function(object,
                                     proper.method = "project",
                                     proper.control = list(),
                                     slice.context = NULL) {
  proper.method <- match.arg(as.character(proper.method)[1L], c("project"))
  proper.control <- .np_condens_normalize_proper_control(proper.control)

  if (!isTRUE(object$trainiseval) && identical(proper.control$apply, "fitted")) {
    info <- .np_condens_make_reason_info(
      reason = "scope_not_selected",
      supported = TRUE,
      slice.count = 0L,
      grid.common = FALSE
    )
    return(list(applied = FALSE, reason = "scope_not_selected", proper.info = info))
  }

  grid.out <- .np_condens_apply_proper_grid(
    object = object,
    proper.method = proper.method,
    proper.control = proper.control
  )
  if (isTRUE(grid.out$applied))
    return(grid.out)

  if (!identical(proper.control$mode, "slice"))
    return(grid.out)

  .np_condens_apply_proper_slice(
    object = object,
    proper.method = proper.method,
    proper.control = proper.control,
    slice.context = slice.context,
    grid.out = grid.out
  )
}

.np_condens_finalize_proper_object <- function(object,
                                               proper = FALSE,
                                               proper.method = c("project"),
                                               proper.control = list(),
                                               slice.context = NULL,
                                               where = "npcdens()") {
  args <- .np_condens_validate_proper_args(
    proper = proper,
    proper.method = proper.method,
    proper.control = proper.control
  )

  object$proper.requested <- args$proper.requested
  object$proper.applied <- FALSE
  object$proper.method <- if (isTRUE(args$proper.requested)) args$proper.method else NULL
  object$condens.raw <- NULL
  object$proper.info <- .np_condens_make_reason_info(
    reason = if (isTRUE(args$proper.requested)) "pending" else "not_requested",
    supported = FALSE
  )

  if (!isTRUE(args$proper.requested))
    return(object)

  if (.np_condens_is_already_proper_by_design(object$bws)) {
    object$proper.info <- .np_condens_make_reason_info(
      reason = "already_proper",
      supported = TRUE
    )
    return(object)
  }

  if (.np_condens_proper_shadow_enabled()) {
    object$proper.info <- .np_condens_make_reason_info(
      reason = "repair_disabled_shadow_mode",
      supported = FALSE
    )
    return(object)
  }

  proper.out <- .np_condens_apply_proper(
    object = object,
    proper.method = args$proper.method,
    proper.control = args$proper.control,
    slice.context = slice.context
  )

  if (!isTRUE(proper.out$applied)) {
    object$proper.info <- proper.out$proper.info
    if (isTRUE(args$proper.control$fail.on.unsupported) &&
        !isTRUE(proper.out$proper.info$supported)) {
      stop(.np_condens_proper_reason_message(
        reason = proper.out$reason,
        where = where
      ), call. = FALSE)
    }
    return(object)
  }

  object$condens <- proper.out$condens
  object$proper.applied <- TRUE
  object$proper.info <- proper.out$proper.info
  if (!is.null(proper.out$condens.raw))
    object$condens.raw <- proper.out$condens.raw

  object
}
