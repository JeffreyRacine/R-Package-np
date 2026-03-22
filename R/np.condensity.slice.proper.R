.np_condens_slice_groups <- function(xeval) {
  xeval <- toFrame(xeval)

  if (ncol(xeval) == 0L) {
    split(seq_len(nrow(xeval)), factor(rep.int("all", nrow(xeval))), drop = TRUE)
  } else {
    groups <- do.call(interaction, c(unname(xeval), list(drop = TRUE, lex.order = TRUE)))
    split(seq_len(nrow(xeval)), groups, drop = TRUE)
  }
}

.np_condens_slice_bounds <- function(object, slice.context, proper.control) {
  yidx <- object$bws$iycon
  if (length(yidx) == 1L &&
      length(object$bws$cykerlb) >= yidx &&
      length(object$bws$cykerub) >= yidx &&
      is.finite(object$bws$cykerlb[yidx]) &&
      is.finite(object$bws$cykerub[yidx])) {
    bounds <- c(object$bws$cykerlb[yidx], object$bws$cykerub[yidx])
  } else {
    bounds <- grDevices::extendrange(
      c(as.double(slice.context$tydat[[1L]]), as.double(slice.context$eydat[[1L]])),
      f = proper.control$`slice.extend.factor`
    )
  }

  bounds <- as.double(bounds)
  if (length(bounds) != 2L ||
      any(!is.finite(bounds)) ||
      !(bounds[2L] > bounds[1L])) {
    stop("invalid slice master-grid bounds")
  }

  bounds
}

.np_condens_build_slice_master_grid <- function(object, slice.context, proper.control) {
  proper.control <- .np_condens_normalize_proper_control(proper.control)

  y.req <- as.double(slice.context$eydat[[1L]])
  if (!length(y.req) || any(!is.finite(y.req)))
    stop("requested evaluation y values must be finite for slice repair")

  bounds <- .np_condens_slice_bounds(
    object = object,
    slice.context = slice.context,
    proper.control = proper.control
  )

  base.grid <- seq(
    from = bounds[1L],
    to = bounds[2L],
    length.out = proper.control$`slice.grid.size`
  )
  y.grid <- sort(unique(c(base.grid, y.req)))

  if (length(y.grid) < 2L || any(diff(y.grid) <= proper.control$tol))
    stop("slice repair master-grid must contain at least two strictly increasing points")

  y.grid
}

.np_condens_build_slice_eval_grid <- function(object, slice.context, proper.control) {
  xeval <- toFrame(slice.context$exdat)
  eydat <- toFrame(slice.context$eydat)
  groups <- .np_condens_slice_groups(xeval)
  y.grid <- .np_condens_build_slice_master_grid(
    object = object,
    slice.context = slice.context,
    proper.control = proper.control
  )

  x.unique <- xeval[vapply(groups, `[`, integer(1), 1L), , drop = FALSE]
  ny <- length(y.grid)
  ng <- length(groups)

  if (ncol(x.unique) == 0L) {
    ex.grid <- data.frame(row.names = seq_len(ng * ny))
  } else {
    ex.grid <- x.unique[rep(seq_len(ng), each = ny), , drop = FALSE]
    rownames(ex.grid) <- NULL
  }

  ey.grid <- eydat[rep.int(1L, ng * ny), , drop = FALSE]
  ey.grid[[1L]] <- rep(y.grid, times = ng)
  rownames(ey.grid) <- NULL

  list(
    groups = groups,
    y.grid = y.grid,
    exdat = ex.grid,
    eydat = ey.grid,
    grid.slices = split(seq_len(ng * ny), rep(seq_len(ng), each = ny))
  )
}

.np_condens_apply_proper_slice <- function(object,
                                           proper.method = "project",
                                           proper.control = list(),
                                           slice.context = NULL,
                                           grid.out = NULL) {
  proper.method <- match.arg(as.character(proper.method)[1L], c("project"))
  proper.control <- .np_condens_normalize_proper_control(proper.control)

  if (!.np_condens_slice_dispatch_enabled()) {
    info <- .np_condens_make_reason_info(
      reason = "slice_disabled",
      supported = FALSE,
      slice.count = 0L,
      grid.common = FALSE
    )
    return(list(applied = FALSE, reason = "slice_disabled", proper.info = info))
  }

  target.context <- NULL
  target.scope <- NULL

  if (isTRUE(object$trainiseval)) {
    if (!identical(proper.control$apply, "fitted") &&
        !identical(proper.control$apply, "both")) {
      return(grid.out)
    }

    if (is.null(slice.context) ||
        is.null(slice.context$txdat) ||
        is.null(slice.context$tydat)) {
      info <- .np_condens_make_reason_info(
        reason = "slice_context_missing",
        supported = FALSE,
        slice.count = 0L,
        grid.common = FALSE
      )
      return(list(applied = FALSE, reason = "slice_context_missing", proper.info = info))
    }

    target.context <- list(
      txdat = slice.context$txdat,
      tydat = slice.context$tydat,
      exdat = slice.context$txdat,
      eydat = slice.context$tydat
    )
    target.scope <- "fitted"
  } else {
    if (identical(proper.control$apply, "fitted")) {
      info <- .np_condens_make_reason_info(
        reason = "scope_not_selected",
        supported = TRUE,
        slice.count = 0L,
        grid.common = FALSE
      )
      return(list(applied = FALSE, reason = "scope_not_selected", proper.info = info))
    }

    if (is.null(slice.context) ||
        is.null(slice.context$txdat) ||
        is.null(slice.context$tydat) ||
        is.null(slice.context$exdat) ||
        is.null(slice.context$eydat)) {
      info <- .np_condens_make_reason_info(
        reason = "slice_context_missing",
        supported = FALSE,
        slice.count = 0L,
        grid.common = FALSE
      )
      return(list(applied = FALSE, reason = "slice_context_missing", proper.info = info))
    }

    target.context <- slice.context
    target.scope <- "evaluation"
  }

  if (is.null(target.context) ||
      is.null(target.context$txdat) ||
      is.null(target.context$tydat) ||
      is.null(target.context$exdat) ||
      is.null(target.context$eydat)) {
    info <- .np_condens_make_reason_info(
      reason = "slice_context_missing",
      supported = FALSE,
      slice.count = 0L,
      grid.common = FALSE
    )
    return(list(applied = FALSE, reason = "slice_context_missing", proper.info = info))
  }

  grid.eval <- tryCatch(
    .np_condens_build_slice_eval_grid(
      object = object,
      slice.context = target.context,
      proper.control = proper.control
    ),
    error = function(e) e
  )
  if (inherits(grid.eval, "error")) {
    info <- .np_condens_make_reason_info(
      reason = "slice_invalid_master_grid",
      supported = FALSE,
      slice.count = 0L,
      grid.common = FALSE
    )
    return(list(
      applied = FALSE,
      reason = "slice_invalid_master_grid",
      proper.info = info
    ))
  }

  grid.fit <- tryCatch(
    npcdens(
      bws = object$bws,
      txdat = target.context$txdat,
      tydat = target.context$tydat,
      exdat = grid.eval$exdat,
      eydat = grid.eval$eydat,
      proper = FALSE
    ),
    error = function(e) e
  )
  if (inherits(grid.fit, "error")) {
    info <- .np_condens_make_reason_info(
      reason = "slice_eval_failed",
      supported = FALSE,
      slice.count = length(grid.eval$groups),
      grid.common = FALSE
    )
    return(list(applied = FALSE, reason = "slice_eval_failed", proper.info = info))
  }

  grid.proper <- .np_condens_apply_proper_grid(
    object = grid.fit,
    proper.method = proper.method,
    proper.control = proper.control
  )
  if (!isTRUE(grid.proper$applied))
    return(grid.proper)

  y.req <- as.double(target.context$eydat[[1L]])
  ypos <- match(y.req, grid.eval$y.grid)
  if (anyNA(ypos)) {
    info <- .np_condens_make_reason_info(
      reason = "slice_invalid_master_grid",
      supported = FALSE,
      slice.count = length(grid.eval$groups),
      grid.common = FALSE
    )
    return(list(
      applied = FALSE,
      reason = "slice_invalid_master_grid",
      proper.info = info
    ))
  }

  repaired <- numeric(length(y.req))
  for (i in seq_along(grid.eval$groups)) {
    idx.req <- grid.eval$groups[[i]]
    idx.grid <- grid.eval$grid.slices[[i]]
    repaired[idx.req] <- grid.proper$condens[idx.grid[ypos[idx.req]]]
  }

  info <- grid.proper$proper.info
  info$route <- "slice"
  info$apply.scope <- target.scope
  info$internal.grid.size <- length(grid.eval$y.grid)
  info$request.nobs <- nrow(target.context$eydat)

  list(
    applied = TRUE,
    condens = repaired,
    condens.raw = if (isTRUE(proper.control$store.raw)) as.double(object$condens) else NULL,
    proper.info = info
  )
}
