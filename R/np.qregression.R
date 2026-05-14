npqreg <-
  function(bws, ...){
    args <- list(...)

    if (!missing(bws)){
      if (is.recursive(bws)){
        if (!is.null(bws$formula) && is.null(args$txdat))
          UseMethod("npqreg",bws$formula)
        else if (!is.null(bws$call) && is.null(args$txdat))
          UseMethod("npqreg",bws$call)
        else if (!is.call(bws))
          UseMethod("npqreg",bws)
        else
          UseMethod("npqreg",NULL)
      } else {
        UseMethod("npqreg", NULL)
      }
    } else {
      UseMethod("npqreg", NULL)
    }
  }

.npqreg.fit.control.names <- c("data", "newdata", "exdat", "tau", "gradients", "tol", "small", "itmax")
.npqreg.removed.solver.controls <- c("ftol",
                                     "lbc.dir", "dfc.dir", "cfac.dir", "initc.dir",
                                     "lbd.dir", "hbd.dir", "dfac.dir", "initd.dir")

.npqreg_validate_tau <- function(tau) {
  if (!is.numeric(tau) || !length(tau) || anyNA(tau) ||
      any(!is.finite(tau)) || any(tau <= 0) || any(tau >= 1))
    stop("'tau' must contain numeric values in (0,1)")
  as.double(tau)
}

.npqreg_tau_labels <- function(tau) {
  paste0("tau=", format(tau, trim = TRUE, scientific = FALSE))
}

.npqreg_napredict_eval <- function(omit, x) {
  if (!length(omit))
    return(x)
  if (length(dim(x)) <= 2L)
    return(napredict(omit, x))
  d <- dim(x)
  dn <- dimnames(x)
  if (!is.null(dn))
    dn[[1L]] <- NULL
  out <- array(NA_real_,
               dim = c(d[1L] + length(omit), d[-1L]),
               dimnames = dn)
  keep <- seq_len(dim(out)[1L])[-as.integer(omit)]
  out[keep, , ] <- x
  out
}

.npqreg_validate_newdata_terms <- function(newdata, xnames) {
  nd <- toFrame(newdata)
  missing.names <- setdiff(xnames, names(nd))
  if (length(missing.names))
    stop(sprintf(
      "newdata must contain columns: %s",
      paste(shQuote(xnames), collapse = ", ")
    ), call. = FALSE)
  invisible(TRUE)
}

.npqreg_fit_dots <- function(dots, allow.bandwidth.controls = FALSE) {
  dot.names <- names(dots)
  if (is.null(dot.names))
    dot.names <- rep("", length(dots))

  stale <- intersect(dot.names[nzchar(dot.names)], .npqreg.removed.solver.controls)
  if (length(stale) && !allow.bandwidth.controls) {
    stop(sprintf(
      "'%s' %s no longer accepted by npqreg; the canonical one-dimensional quantile extractor is controlled by 'tol', 'small', and 'itmax'",
      paste(stale, collapse = "', '"),
      if (length(stale) == 1L) "is" else "are"
    ))
  }

  if (!allow.bandwidth.controls) {
    bad <- dot.names == "" | !(dot.names %in% .npqreg.fit.control.names)
    if (any(bad))
      .np_reject_unused_dots(dots[bad], "npqreg")
  }

  keep <- (!nzchar(dot.names)) | (dot.names %in% .npqreg.fit.control.names)
  dots[keep]
}

.npqreg_strip_fit_controls_from_bw_call <- function(call) {
  for (nm in c("tau", "gradients", "tol", "small", "itmax", "newdata", "exdat")) {
    if (nm %in% names(call))
      call[[nm]] <- NULL
  }
  call
}

.npqreg_quantile_delta_from_conditional <- function(bws,
                                                    xdat,
                                                    ydat,
                                                    exdat,
                                                    quantile,
                                                    gradients = FALSE) {
  xdat <- toFrame(xdat)
  ydat <- toFrame(ydat)
  exdat <- toFrame(exdat)
  quantile <- as.double(quantile)
  gradients <- npValidateScalarLogical(gradients, "gradients")

  if (length(quantile) != nrow(exdat))
    stop("quantile delta helper requires one quantile per evaluation row")
  if (ncol(ydat) != 1L)
    stop("quantile delta helper requires a single response")

  eydat <- stats::setNames(data.frame(quantile), names(ydat)[1L])
  cdf.obj <- .np_conditional_eval_selected(
    bws = bws,
    xdat = xdat,
    ydat = ydat,
    exdat = exdat,
    eydat = eydat,
    cdf = TRUE,
    gradients = gradients
  )
  dens.obj <- .np_conditional_eval_selected(
    bws = bws,
    xdat = xdat,
    ydat = ydat,
    exdat = exdat,
    eydat = eydat,
    cdf = FALSE,
    gradients = FALSE
  )

  dens <- as.double(dens.obj$condens)
  quanterr <- as.double(cdf.obj$conderr) / NZD(dens)
  quanterr[!is.finite(quanterr) | quanterr < 0.0] <- NA_real_

  if (!gradients) {
    return(list(
      quanterr = quanterr,
      quantgrad = NA,
      quantgerr = NA,
      cdf = cdf.obj,
      dens = dens.obj
    ))
  }

  dens.mat <- matrix(NZD(dens),
                     nrow = nrow(cdf.obj$congrad),
                     ncol = ncol(cdf.obj$congrad))
  grad <- -cdf.obj$congrad / dens.mat
  grad[!is.finite(grad)] <- NA_real_

  gerr <- cdf.obj$congerr / dens.mat
  gerr[!is.finite(gerr) | gerr < 0.0] <- NA_real_

  list(
    quanterr = quanterr,
    quantgrad = grad,
    quantgerr = gerr,
    cdf = cdf.obj,
    dens = dens.obj
  )
}

.npqreg_selected_cdf_values <- function(bws,
                                        xdat,
                                        ydat,
                                        exdat,
                                        ycand) {
  ydat <- toFrame(ydat)
  yname <- names(ydat)[1L]
  eydat <- stats::setNames(data.frame(as.double(ycand)), yname)
  as.double(.np_conditional_eval_selected(
    bws = bws,
    xdat = xdat,
    ydat = ydat,
    exdat = exdat,
    eydat = eydat,
    cdf = TRUE,
    gradients = FALSE
  )$condist)
}

.npRmpi_npqreg_parallel_context <- function(bws, comm = 1L) {
  isTRUE(isa(bws, "condbandwidth")) &&
    isTRUE(.npRmpi_has_active_slave_pool(comm = comm)) &&
    !isTRUE(getOption("npRmpi.local.regression.mode", FALSE)) &&
    !isTRUE(.npRmpi_autodispatch_called_from_bcast())
}

.npRmpi_npqreg_chunk_size <- function(n.eval, comm = 1L) {
  n.eval <- as.integer(n.eval)
  if (is.na(n.eval) || n.eval < 1L)
    return(1L)

  opt <- suppressWarnings(as.integer(getOption("npRmpi.npqreg.chunk.size", NA_integer_))[1L])
  if (!is.na(opt) && opt > 0L)
    return(min(n.eval, opt))

  workers <- .npRmpi_bootstrap_worker_count(comm = comm)
  slots <- max(1L, workers + 1L)
  max(1L, as.integer(ceiling(n.eval / slots)))
}

.npRmpi_npqreg_parallel_min_eval <- function(comm = 1L) {
  opt <- suppressWarnings(as.integer(getOption("npRmpi.npqreg.parallel.min.eval", NA_integer_))[1L])
  if (!is.na(opt) && opt > 0L)
    return(opt)

  workers <- .npRmpi_bootstrap_worker_count(comm = comm)
  max(2000L, 8L * (workers + 1L))
}

.npRmpi_npqreg_parallel_ready <- function(bws,
                                          n.eval,
                                          comm = 1L,
                                          what = "npqreg",
                                          force.parallel = FALSE) {
  n.eval <- as.integer(n.eval)
  if (is.na(n.eval))
    return(FALSE)
  if (!isTRUE(force.parallel) &&
      n.eval < .npRmpi_npqreg_parallel_min_eval(comm = comm))
    return(FALSE)
  if (!.npRmpi_npqreg_parallel_context(bws, comm = comm))
    return(FALSE)

  .npRmpi_bootstrap_fanout_enabled(
    comm = comm,
    n = n.eval,
    B = n.eval,
    chunk.size = .npRmpi_npqreg_chunk_size(n.eval = n.eval, comm = comm),
    what = what
  )
  TRUE
}

.npRmpi_npqreg_reset_worker_comm_state <- function(comm = 1L) {
  if (isTRUE(.npRmpi_autodispatch_called_from_bcast()) ||
      !isTRUE(.npRmpi_has_active_slave_pool(comm = comm)))
    return(invisible(FALSE))

  cmd <- quote({
    try(.Call("C_np_set_local_regression_mode", FALSE, PACKAGE = "npRmpi"),
        silent = TRUE)
    try(.Call("C_np_set_active_comm", FALSE, as.integer(1L), PACKAGE = "npRmpi"),
        silent = TRUE)
    invisible(NULL)
  })
  try(.npRmpi_bcast_cmd_expr(cmd, comm = comm, caller.execute = TRUE),
      silent = TRUE)
  invisible(TRUE)
}

.npqreg_selected_cdf_values_parallel <- function(bws,
                                                 xdat,
                                                 ydat,
                                                 exdat,
                                                 ycand,
                                                 comm = 1L) {
  exdat <- toFrame(exdat)
  n.eval <- nrow(exdat)
  if (!.npRmpi_npqreg_parallel_ready(
        bws = bws,
        n.eval = n.eval,
        comm = comm,
        what = "npqreg selected CDF"
      )) {
    return(.npRmpi_with_local_cdist_eval(.npqreg_selected_cdf_values(
      bws = bws,
      xdat = xdat,
      ydat = ydat,
      exdat = exdat,
      ycand = ycand
    )))
  }

  ycand <- as.double(ycand)
  on.exit(.npRmpi_npqreg_reset_worker_comm_state(comm = comm), add = TRUE)
  tasks <- .npRmpi_bootstrap_chunk_tasks(
    B = n.eval,
    chunk.size = .npRmpi_npqreg_chunk_size(n.eval = n.eval, comm = comm)
  )
  worker <- function(task, bws, xdat, ydat, exdat, ycand) {
    idx <- seq.int(as.integer(task$start),
                   length.out = as.integer(task$bsz))
    .npRmpi_with_local_cdist_eval(.npqreg_selected_cdf_values(
      bws = bws,
      xdat = xdat,
      ydat = ydat,
      exdat = exdat[idx, , drop = FALSE],
      ycand = ycand[idx]
    ))
  }

  out <- .npRmpi_bootstrap_run_fanout(
    tasks = tasks,
    worker = worker,
    ncol.out = 1L,
    what = "npqreg selected CDF",
    progress.label = "npqreg selected CDF",
    profile.where = "npqreg:selected-cdf",
    comm = comm,
    master_local_chunk = TRUE,
    bws = bws,
    xdat = xdat,
    ydat = ydat,
    exdat = exdat,
    ycand = ycand
  )

  as.double(out[, 1L])
}

.npqreg_quantile_delta_matrix <- function(delta, gradients = FALSE) {
  quanterr <- as.double(delta$quanterr)
  if (!isTRUE(gradients))
    return(matrix(quanterr, ncol = 1L))

  cbind(
    quanterr,
    as.matrix(delta$quantgrad),
    as.matrix(delta$quantgerr)
  )
}

.npqreg_quantile_delta_from_conditional_parallel <- function(bws,
                                                             xdat,
                                                             ydat,
                                                             exdat,
                                                             quantile,
                                                             gradients = FALSE,
                                                             comm = 1L) {
  exdat <- toFrame(exdat)
  n.eval <- nrow(exdat)
  gradients <- npValidateScalarLogical(gradients, "gradients")
  grad.cols <- if (isTRUE(gradients)) as.integer(bws$xndim) else 0L
  if (is.na(grad.cols) || grad.cols < 0L)
    grad.cols <- 0L

  if (!.npRmpi_npqreg_parallel_ready(
        bws = bws,
        n.eval = n.eval,
        comm = comm,
        what = "npqreg quantile delta"
      )) {
    return(.npRmpi_with_local_cdist_eval(.npqreg_quantile_delta_from_conditional(
      bws = bws,
      xdat = xdat,
      ydat = ydat,
      exdat = exdat,
      quantile = quantile,
      gradients = gradients
    )))
  }

  quantile <- as.double(quantile)
  on.exit(.npRmpi_npqreg_reset_worker_comm_state(comm = comm), add = TRUE)
  tasks <- .npRmpi_bootstrap_chunk_tasks(
    B = n.eval,
    chunk.size = .npRmpi_npqreg_chunk_size(n.eval = n.eval, comm = comm)
  )
  worker <- function(task, bws, xdat, ydat, exdat, quantile, gradients) {
    idx <- seq.int(as.integer(task$start),
                   length.out = as.integer(task$bsz))
    .npRmpi_with_local_cdist_eval(.npqreg_quantile_delta_matrix(
      .npqreg_quantile_delta_from_conditional(
        bws = bws,
        xdat = xdat,
        ydat = ydat,
        exdat = exdat[idx, , drop = FALSE],
        quantile = quantile[idx],
        gradients = gradients
      ),
      gradients = gradients
    ))
  }

  out <- .npRmpi_bootstrap_run_fanout(
    tasks = tasks,
    worker = worker,
    ncol.out = 1L + 2L * grad.cols,
    what = "npqreg quantile delta",
    progress.label = "npqreg quantile delta",
    profile.where = "npqreg:quantile-delta",
    comm = comm,
    master_local_chunk = TRUE,
    bws = bws,
    xdat = xdat,
    ydat = ydat,
    exdat = exdat,
    quantile = quantile,
    gradients = gradients
  )

  quanterr <- as.double(out[, 1L])
  if (!isTRUE(gradients)) {
    return(list(
      quanterr = quanterr,
      quantgrad = NA,
      quantgerr = NA,
      cdf = NULL,
      dens = NULL
    ))
  }

  grad.idx <- seq.int(2L, length.out = grad.cols)
  gerr.idx <- seq.int(2L + grad.cols, length.out = grad.cols)
  list(
    quanterr = quanterr,
    quantgrad = out[, grad.idx, drop = FALSE],
    quantgerr = out[, gerr.idx, drop = FALSE],
    cdf = NULL,
    dens = NULL
  )
}

.npqreg_fit_tau_vector_parallel_matrix <- function(bws,
                                                   xdat,
                                                   ydat,
                                                   exdat,
                                                   tau,
                                                   gradients = FALSE,
                                                   tol,
                                                   small,
                                                   itmax,
                                                   comm = 1L,
                                                   force.parallel = FALSE) {
  exdat <- toFrame(exdat)
  n.eval <- nrow(exdat)
  tau <- .npqreg_validate_tau(tau)
  gradients <- npValidateScalarLogical(gradients, "gradients")
  grad.cols <- if (isTRUE(gradients)) as.integer(bws$xndim) else 0L
  if (is.na(grad.cols) || grad.cols < 0L)
    grad.cols <- 0L
  cols.per.tau <- 2L + 2L * grad.cols

  fit_chunk <- function(ex.chunk) {
    pieces <- vector("list", length(tau))
    for (j in seq_along(tau)) {
      yq <- .npRmpi_with_local_cdist_eval(.npqreg_invert_selected_cdf(
        bws = bws,
        xdat = xdat,
        ydat = ydat,
        exdat = ex.chunk,
        tau = tau[[j]],
        tol = tol,
        small = small,
        itmax = itmax,
        parallel = FALSE
      ))
      delta <- .npRmpi_with_local_cdist_eval(.npqreg_quantile_delta_from_conditional(
        bws = bws,
        xdat = xdat,
        ydat = ydat,
        exdat = ex.chunk,
        quantile = yq,
        gradients = gradients
      ))
      pieces[[j]] <- cbind(yq, .npqreg_quantile_delta_matrix(delta, gradients = gradients))
    }
    do.call(cbind, pieces)
  }

  if (!.npRmpi_npqreg_parallel_ready(
        bws = bws,
        n.eval = n.eval,
        comm = comm,
        what = "npqreg tau block",
        force.parallel = force.parallel
      )) {
    return(fit_chunk(exdat))
  }

  on.exit(.npRmpi_npqreg_reset_worker_comm_state(comm = comm), add = TRUE)
  tasks <- .npRmpi_bootstrap_chunk_tasks(
    B = n.eval,
    chunk.size = .npRmpi_npqreg_chunk_size(n.eval = n.eval, comm = comm)
  )
  worker <- function(task, bws, xdat, ydat, exdat, tau, gradients, tol, small, itmax) {
    idx <- seq.int(as.integer(task$start),
                   length.out = as.integer(task$bsz))
    ex.chunk <- exdat[idx, , drop = FALSE]
    pieces <- vector("list", length(tau))
    for (j in seq_along(tau)) {
      yq <- .npRmpi_with_local_cdist_eval(.npqreg_invert_selected_cdf(
        bws = bws,
        xdat = xdat,
        ydat = ydat,
        exdat = ex.chunk,
        tau = tau[[j]],
        tol = tol,
        small = small,
        itmax = itmax,
        parallel = FALSE
      ))
      delta <- .npRmpi_with_local_cdist_eval(.npqreg_quantile_delta_from_conditional(
        bws = bws,
        xdat = xdat,
        ydat = ydat,
        exdat = ex.chunk,
        quantile = yq,
        gradients = gradients
      ))
      pieces[[j]] <- cbind(yq, .npqreg_quantile_delta_matrix(delta, gradients = gradients))
    }
    do.call(cbind, pieces)
  }

  .npRmpi_bootstrap_run_fanout(
    tasks = tasks,
    worker = worker,
    ncol.out = length(tau) * cols.per.tau,
    what = "npqreg tau block",
    progress.label = "npqreg tau block",
    profile.where = "npqreg:tau-block",
    comm = comm,
    master_local_chunk = TRUE,
    bws = bws,
    xdat = xdat,
    ydat = ydat,
    exdat = exdat,
    tau = tau,
    gradients = gradients,
    tol = tol,
    small = small,
    itmax = itmax
  )
}

.npqreg_fit_tau_vector_from_parallel_matrix <- function(mat,
                                                        tau,
                                                        gradients = FALSE,
                                                        grad.names = NULL) {
  tau <- .npqreg_validate_tau(tau)
  gradients <- npValidateScalarLogical(gradients, "gradients")
  grad.cols <- if (isTRUE(gradients)) {
    per.tau.raw <- ncol(mat) / length(tau)
    (per.tau.raw - 2L) / 2L
  } else {
    0L
  }
  if (!is.finite(grad.cols) || grad.cols < 0L || grad.cols != floor(grad.cols))
    stop("internal error: malformed npqreg parallel gradient payload", call. = FALSE)
  grad.cols <- as.integer(grad.cols)
  cols.per.tau <- 2L + 2L * grad.cols
  if (ncol(mat) != length(tau) * cols.per.tau)
    stop("internal error: malformed npqreg parallel tau payload", call. = FALSE)

  n.eval <- nrow(mat)
  yq <- matrix(NA_real_, nrow = n.eval, ncol = length(tau))
  yqerr <- matrix(NA_real_, nrow = n.eval, ncol = length(tau))
  if (isTRUE(gradients)) {
    yqgrad <- array(NA_real_, dim = c(n.eval, grad.cols, length(tau)))
    yqgerr <- array(NA_real_, dim = c(n.eval, grad.cols, length(tau)))
  }

  for (j in seq_along(tau)) {
    offset <- (j - 1L) * cols.per.tau
    yq[, j] <- mat[, offset + 1L]
    yqerr[, j] <- mat[, offset + 2L]
    if (isTRUE(gradients) && grad.cols > 0L) {
      grad.idx <- seq.int(offset + 3L, length.out = grad.cols)
      gerr.idx <- seq.int(offset + 3L + grad.cols, length.out = grad.cols)
      yqgrad[, , j] <- mat[, grad.idx, drop = FALSE]
      yqgerr[, , j] <- mat[, gerr.idx, drop = FALSE]
    }
  }

  tau.labels <- .npqreg_tau_labels(tau)
  if (length(tau) == 1L) {
    return(list(
      yq = as.double(yq[, 1L]),
      yqerr = as.double(yqerr[, 1L]),
      yqgrad = if (isTRUE(gradients)) {
        out <- yqgrad[, , 1L, drop = FALSE]
        dim(out) <- c(n.eval, grad.cols)
        if (!is.null(grad.names) && length(grad.names) == grad.cols)
          colnames(out) <- grad.names
        out
      } else NA,
      yqgerr = if (isTRUE(gradients)) {
        out <- yqgerr[, , 1L, drop = FALSE]
        dim(out) <- c(n.eval, grad.cols)
        if (!is.null(grad.names) && length(grad.names) == grad.cols)
          colnames(out) <- grad.names
        out
      } else NA
    ))
  }

  colnames(yq) <- tau.labels
  colnames(yqerr) <- tau.labels
  if (isTRUE(gradients)) {
    dimnames(yqgrad) <- list(NULL, NULL, tau.labels)
    dimnames(yqgerr) <- list(NULL, NULL, tau.labels)
    if (!is.null(grad.names) && length(grad.names) == grad.cols) {
      dimnames(yqgrad)[[2L]] <- grad.names
      dimnames(yqgerr)[[2L]] <- grad.names
    }
  }
  list(
    yq = yq,
    yqerr = yqerr,
    yqgrad = if (isTRUE(gradients)) yqgrad else NA,
    yqgerr = if (isTRUE(gradients)) yqgerr else NA
  )
}

.npqreg_assert_selected_cdf_metadata <- function(bws) {
  reg.engine <- if (is.null(bws$regtype.engine)) {
    if (is.null(bws$regtype)) "lc" else as.character(bws$regtype)
  } else {
    as.character(bws$regtype.engine)
  }
  if (!identical(reg.engine, "lp"))
    return(invisible(TRUE))

  if (is.null(bws$degree.engine) && is.null(bws$degree))
    stop("selected LP conditional distribution metadata missing from bandwidth object: degree")
  invisible(TRUE)
}

.npqreg_invert_selected_cdf <- function(bws,
                                        xdat,
                                        ydat,
                                        exdat,
                                        tau,
                                        tol,
                                        small,
                                        itmax,
                                        parallel = FALSE,
                                        comm = 1L) {
  .npqreg_assert_selected_cdf_metadata(bws)

  xdat <- toFrame(xdat)
  ydat <- toFrame(ydat)
  exdat <- toFrame(exdat)
  y <- as.double(ydat[[1L]])
  y <- y[is.finite(y)]
  if (!length(y))
    stop("npqreg selected-CDF inversion requires finite response support")

  n.eval <- nrow(exdat)
  y.min <- min(y)
  y.max <- max(y)
  if (!is.finite(y.min) || !is.finite(y.max))
    stop("npqreg selected-CDF inversion found non-finite response support")
  if (identical(y.min, y.max))
    return(rep.int(y.min, n.eval))

  lo <- rep.int(y.min, n.eval)
  hi <- rep.int(y.max, n.eval)
  cdf_values <- if (isTRUE(parallel)) {
    function(bws, xdat, ydat, exdat, ycand) {
      .npqreg_selected_cdf_values_parallel(
        bws = bws,
        xdat = xdat,
        ydat = ydat,
        exdat = exdat,
        ycand = ycand,
        comm = comm
      )
    }
  } else {
    function(bws, xdat, ydat, exdat, ycand) {
      .npRmpi_with_local_cdist_eval(.npqreg_selected_cdf_values(
        bws = bws,
        xdat = xdat,
        ydat = ydat,
        exdat = exdat,
        ycand = ycand
      ))
    }
  }

  flo <- cdf_values(bws, xdat, ydat, exdat, lo)
  fhi <- cdf_values(bws, xdat, ydat, exdat, hi)
  if (any(!is.finite(flo)) || any(!is.finite(fhi)))
    stop("npqreg selected-CDF inversion encountered non-finite bracket values")

  done.low <- flo >= tau
  done.high <- fhi < tau
  active <- !(done.low | done.high)

  maxiter <- min(as.integer(itmax), 1000L)
  iter <- 0L
  while (any(active) && iter < maxiter) {
    iter <- iter + 1L
    mid <- (lo[active] + hi[active]) / 2.0
    fmid <- cdf_values(
      bws = bws,
      xdat = xdat,
      ydat = ydat,
      exdat = exdat[active, , drop = FALSE],
      ycand = mid
    )
    if (any(!is.finite(fmid)))
      stop("npqreg selected-CDF inversion encountered non-finite refinement values")

    active.idx <- which(active)
    upper <- fmid >= tau
    hi[active.idx[upper]] <- mid[upper]
    lo[active.idx[!upper]] <- mid[!upper]

    width <- hi[active.idx] - lo[active.idx]
    scale <- pmax(abs(hi[active.idx]), abs(lo[active.idx]), 1.0)
    active[active.idx] <- width > (tol * scale + small)
  }

  if (any(active))
    stop("npqreg selected-CDF inversion failed to converge within 'itmax'")

  out <- hi
  out[done.low] <- y.min
  out[done.high] <- y.max
  out
}

npqreg.formula <-
  function(bws, data = NULL, newdata = NULL, ...){

    tt <- terms(bws)
    m <- match(c("formula", "data", "subset", "na.action"),
               names(bws$call), nomatch = 0)
    tmf <- bws$call[c(1,m)]
    tmf[[1]] <- as.name("model.frame")
    tmf[["formula"]] <- tt
    mf.args <- as.list(tmf)[-1L]
    umf <- tmf <- do.call(stats::model.frame, mf.args, envir = environment(tt))

    tydat <- tmf[, bws$variableNames[["response"]], drop = FALSE]
    txdat <- tmf[, bws$variableNames[["terms"]], drop = FALSE]

    has.eval <- !is.null(newdata)
    if (has.eval) {
      .npqreg_validate_newdata_terms(newdata, bws$variableNames[["terms"]])
      tt <- drop.terms(tt, match(bws$variableNames$response, attr(tt, 'term.labels')))
      umf.args <- list(formula = tt, data = newdata)
      umf <- do.call(stats::model.frame, umf.args, envir = parent.frame())
      emf <- umf
      exdat <- emf[, bws$variableNames[["terms"]], drop = FALSE]
    }

    q.args <- list(txdat = txdat, tydat = tydat)
    if (has.eval)
      q.args$exdat <- exdat
    q.args$bws <- bws
    tbw <- do.call(npqreg, c(q.args, .npqreg_fit_dots(list(...))))

    tbw$omit <- attr(umf,"na.action")
    tbw$rows.omit <- as.vector(tbw$omit)
    tbw$nobs.omit <- length(tbw$rows.omit)

    tbw$quantile <- .npqreg_napredict_eval(tbw$omit, tbw$quantile)
    tbw$quanterr <- .npqreg_napredict_eval(tbw$omit, tbw$quanterr)

    if(tbw$gradients){
        tbw$quantgrad <- .npqreg_napredict_eval(tbw$omit, tbw$quantgrad)
        tbw$quantgerr <- .npqreg_napredict_eval(tbw$omit, tbw$quantgerr)
    }

    return(tbw)
  }

npqreg.call <-
  function(bws, ...) {
    npqreg(txdat = .np_eval_bws_call_arg(bws, "xdat"),
           tydat = .np_eval_bws_call_arg(bws, "ydat"),
           bws = bws, ...)
  }

npqreg.conbandwidth <-
  function(bws, ...){
    stop("incorrect bandwidth type: expected conditional distribution bandwidths instead of conditional density bandwidths")
  }

.npRmpi_npqreg_should_localize <- function(bws) {
  isa(bws, "condbandwidth")
}

.npRmpi_npqreg_eval_local_no_dispatch <- function(expr) {
  old.disable <- getOption("npRmpi.autodispatch.disable", FALSE)
  old.ctx <- getOption("npRmpi.autodispatch.context", FALSE)
  options(npRmpi.autodispatch.disable = TRUE)
  options(npRmpi.autodispatch.context = TRUE)
  on.exit(options(npRmpi.autodispatch.disable = old.disable), add = TRUE)
  on.exit(options(npRmpi.autodispatch.context = old.ctx), add = TRUE)
  on.exit(.npRmpi_npqreg_reset_worker_comm_state(comm = 1L), add = TRUE)
  force(expr)
}

npqreg.condbandwidth <-
  function(bws,
           txdat = stop("training data 'txdat' missing"),
           tydat = stop("training data 'tydat' missing"),
           exdat,
           tau = 0.5,
           gradients = FALSE,
           tol = 1.490116e-04,
           small = 1.490116e-05, itmax = 10000,
           ...){

    fit.start <- proc.time()[3]
    tau <- .npqreg_validate_tau(tau)
    fit.dots <- .npqreg_fit_dots(list(...))
    if (length(fit.dots))
      stop(sprintf("unused npqreg fit argument '%s'", names(fit.dots)[1L]))
    gradients <- npValidateScalarLogical(gradients, "gradients")
    if (!is.numeric(itmax) || length(itmax) != 1L || is.na(itmax) ||
        !is.finite(itmax) || itmax < 1 || itmax != floor(itmax))
      stop("'itmax' must be a positive integer")
    if (!is.numeric(tol) || length(tol) != 1L || is.na(tol) ||
        !is.finite(tol) || tol <= 0)
      stop("'tol' must be a positive finite numeric scalar")
    if (!is.numeric(small) || length(small) != 1L || is.na(small) ||
        !is.finite(small) || small <= 0)
      stop("'small' must be a positive finite numeric scalar")
    itmax <- as.integer(itmax)
    tol <- as.double(tol)
    small <- as.double(small)
    .npRmpi_require_active_slave_pool(where = "npqreg()")
    parallel.cond <- .npRmpi_npqreg_parallel_context(bws, comm = 1L)
    if (isTRUE(parallel.cond)) {
      n.eval.pre <- suppressWarnings(as.integer(tryCatch(
        if (missing(exdat)) NROW(txdat) else NROW(exdat),
        error = function(e) NA_integer_
      ))[1L])
      if (is.na(n.eval.pre) ||
          n.eval.pre < .npRmpi_npqreg_parallel_min_eval(comm = 1L))
        parallel.cond <- FALSE
    }
    if (isTRUE(parallel.cond))
      on.exit(.npRmpi_npqreg_reset_worker_comm_state(comm = 1L), add = TRUE)
    if (.npRmpi_npqreg_should_localize(bws) &&
        !isTRUE(getOption("npRmpi.local.regression.mode", FALSE)) &&
        !isTRUE(.npRmpi_autodispatch_in_context()) &&
        !isTRUE(parallel.cond))
      return(.npRmpi_npqreg_eval_local_no_dispatch(
        .npRmpi_eval_without_dispatch(match.call(), parent.frame())
      ))
    if (.npRmpi_autodispatch_active() && !isTRUE(parallel.cond))
      return(.npRmpi_autodispatch_call(match.call(), parent.frame()))

    no.ex = missing(exdat)

    txdat = toFrame(txdat)
    tydat = toFrame(tydat)

    ntau <- length(tau)
    tau.labels <- .npqreg_tau_labels(tau)

    if (dim(tydat)[2] != 1)
      stop("'tydat' has more than one column")

    if (!no.ex){
      exdat = toFrame(exdat)
      
      if (! txdat %~% exdat )
        stop("'txdat' and 'exdat' are not similar data frames!")
    }

    if (length(bws$xbw) != length(txdat))
      stop("length of bandwidth vector does not match number of columns of 'txdat'")

    if (length(bws$ybw) != 1)
      stop("length of bandwidth vector does not match number of columns of 'tydat'")

    if (any(bws$iyord) || any(bws$iyuno) || coarseclass(tydat[,1]) != "numeric")
      stop("'tydat' is not continuous")

    if ((any(bws$ixcon) &&
         !all(vapply(txdat[, bws$ixcon, drop = FALSE], inherits, logical(1), c("integer", "numeric")))) ||
        (any(bws$ixord) &&
         !all(vapply(txdat[, bws$ixord, drop = FALSE], inherits, logical(1), "ordered"))) ||
        (any(bws$ixuno) &&
         !all(vapply(txdat[, bws$ixuno, drop = FALSE], inherits, logical(1), "factor"))))
      stop("supplied bandwidths do not match 'txdat' in type")

    ## catch and destroy NA's
    keep.rows <- rep_len(TRUE, nrow(txdat))
    rows.omit <- attr(na.omit(data.frame(txdat, tydat)), "na.action")
    if (length(rows.omit) > 0L)
      keep.rows[as.integer(rows.omit)] <- FALSE

    if (!any(keep.rows))
      stop("Training data has no rows without NAs")

    txdat <- txdat[keep.rows,,drop = FALSE]
    tydat <- tydat[keep.rows,,drop = FALSE]

    if (!no.ex){
      keep.eval <- rep_len(TRUE, nrow(exdat))
      rows.omit <- attr(na.omit(exdat), "na.action")
      if (length(rows.omit) > 0L)
        keep.eval[as.integer(rows.omit)] <- FALSE
      exdat <- exdat[keep.eval,,drop = FALSE]
    }
    
    tnrow = dim(txdat)[1]
    enrow = (if (no.ex) tnrow else dim(exdat)[1])

    ## re-assign levels in training and evaluation data to ensure correct
    ## conversion to numeric type.
    
    txdat <- adjustLevels(txdat, bws$xdati)
    tydat <- adjustLevels(tydat, bws$ydati)
    
    if (!no.ex){
      exdat <- adjustLevels(exdat, bws$xdati)
    }

    ## grab the evaluation data before it is converted to numeric
    if(no.ex){
      txeval <- txdat
    } else {
      txeval <- exdat
    }
    txdat.df <- txdat
    tydat.df <- tydat
    if (!no.ex)
      exdat.df <- exdat

    if (isTRUE(parallel.cond)) {
      mat <- .npqreg_fit_tau_vector_parallel_matrix(
        bws = bws,
        xdat = txdat.df,
        ydat = tydat.df,
        exdat = txeval,
        tau = tau,
        gradients = gradients,
        tol = tol,
        small = small,
        itmax = itmax,
        comm = 1L,
        force.parallel = TRUE
      )
      myout <- .npqreg_fit_tau_vector_from_parallel_matrix(
        mat,
        tau = tau,
        gradients = gradients
      )
    } else {
      fit_one_tau <- function(tau_i) {
        yq <- .npqreg_invert_selected_cdf(
          bws = bws,
          xdat = txdat.df,
          ydat = tydat.df,
          exdat = txeval,
          tau = tau_i,
          tol = tol,
          small = small,
          itmax = itmax
        )
        qdelta <- .npRmpi_with_local_cdist_eval(
          .npqreg_quantile_delta_from_conditional(
            bws = bws,
            xdat = txdat.df,
            ydat = tydat.df,
            exdat = txeval,
            quantile = yq,
            gradients = gradients
          )
        )
        list(
          yq = yq,
          yqerr = qdelta$quanterr,
          yqgrad = if (gradients) qdelta$quantgrad else NA,
          yqgerr = if (gradients) qdelta$quantgerr else NA
        )
      }

      tau.out <- lapply(tau, fit_one_tau)

      if (ntau == 1L) {
        myout <- tau.out[[1L]]
      } else {
        myout <- list(
          yq = do.call(cbind, lapply(tau.out, `[[`, "yq")),
          yqerr = do.call(cbind, lapply(tau.out, `[[`, "yqerr")),
          yqgrad = NA,
          yqgerr = NA
        )
        colnames(myout$yq) <- tau.labels
        colnames(myout$yqerr) <- tau.labels
        if (gradients) {
          p <- ncol(tau.out[[1L]]$yqgrad)
          grad.names <- colnames(tau.out[[1L]]$yqgrad)
          myout$yqgrad <- array(NA_real_,
                                dim = c(enrow, p, ntau),
                                dimnames = list(NULL, grad.names, tau.labels))
          myout$yqgerr <- array(NA_real_,
                                dim = c(enrow, p, ntau),
                                dimnames = list(NULL, grad.names, tau.labels))
          for (j in seq_len(ntau)) {
            myout$yqgrad[, , j] <- tau.out[[j]]$yqgrad
            myout$yqgerr[, , j] <- tau.out[[j]]$yqgerr
          }
        }
      }
    }


    fit.elapsed <- proc.time()[3] - fit.start
    optim.time <- if (!is.null(bws$total.time) && is.finite(bws$total.time)) as.double(bws$total.time) else NA_real_
    total.time <- fit.elapsed + (if (is.na(optim.time)) 0.0 else optim.time)

    qregression(bws = bws,
                xeval = txeval,
                tau = tau,
                quantile = myout$yq,
                quanterr = myout$yqerr,
                quantgrad = myout$yqgrad,
                quantgerr = myout$yqgerr,
                ntrain = tnrow,
                trainiseval = no.ex,
                gradients = gradients,
                timing = bws$timing, total.time = total.time,
                optim.time = optim.time, fit.time = fit.elapsed)
  }


npqreg.default <- function(bws, txdat, tydat, nomad = FALSE, ...){
  nomad <- npValidateScalarLogical(nomad, "nomad")
  early.dots <- list(...)
  if ("tau" %in% names(early.dots))
    .npqreg_validate_tau(early.dots$tau)

  if (!missing(bws) && inherits(bws, "formula")) {
    dots <- list(...)
    dot.names <- names(dots)
    if (is.null(dot.names))
      dot.names <- rep("", length(dots))
    fit.names <- c("newdata", "exdat", "tau", "gradients", "tol", "small", "itmax")
    fit.dots <- .npqreg_fit_dots(dots[nzchar(dot.names) & dot.names %in% fit.names])
    bw.dots <- dots[!(nzchar(dot.names) & dot.names %in% fit.names)]
    bw.args <- c(list(formula = bws, nomad = nomad), bw.dots)
    bw.call <- as.call(c(list(quote(npcdistbw)), bw.args))
    use.outer.bandwidth.progress <- !.np_bw_call_uses_nomad_degree_search(
      bw.call,
      caller_env = parent.frame()
    )
    tbw <- if (use.outer.bandwidth.progress) {
      .np_progress_select_bandwidth_enhanced(
        "Selecting conditional distribution bandwidth",
        do.call(npcdistbw, bw.args)
      )
    } else {
      do.call(npcdistbw, bw.args)
    }
    return(do.call(npqreg, c(list(bws = tbw), fit.dots)))
  }

  if (!missing(txdat) && inherits(txdat, "formula") &&
      !missing(bws) && !isa(bws, "condbandwidth")) {
    dots <- list(...)
    dot.names <- names(dots)
    if (is.null(dot.names))
      dot.names <- rep("", length(dots))
    fit.names <- c("newdata", "exdat", "tau", "gradients", "tol", "small", "itmax")
    fit.dots <- .npqreg_fit_dots(dots[nzchar(dot.names) & dot.names %in% fit.names])
    bw.dots <- dots[!(nzchar(dot.names) & dot.names %in% fit.names)]
    bw.args <- c(list(formula = txdat, bws = bws, nomad = nomad), bw.dots)
    bw.call <- as.call(c(list(quote(npcdistbw)), bw.args))
    use.outer.bandwidth.progress <- !.np_bw_call_uses_nomad_degree_search(
      bw.call,
      caller_env = parent.frame()
    )
    tbw <- if (use.outer.bandwidth.progress) {
      .np_progress_select_bandwidth_enhanced(
        "Selecting conditional distribution bandwidth",
        do.call(npcdistbw, bw.args)
      )
    } else {
      do.call(npcdistbw, bw.args)
    }
    return(do.call(npqreg, c(list(bws = tbw), fit.dots)))
  }

  .npRmpi_require_active_slave_pool(where = "npqreg()")
  parallel.cond <- (!missing(bws)) &&
    .npRmpi_npqreg_should_localize(bws) &&
    .npRmpi_npqreg_parallel_context(bws, comm = 1L)
  if (isTRUE(parallel.cond)) {
    n.eval.pre <- suppressWarnings(as.integer(tryCatch(
      if (missing(txdat)) NA_integer_ else NROW(txdat),
      error = function(e) NA_integer_
    ))[1L])
    if (is.na(n.eval.pre) ||
        n.eval.pre < .npRmpi_npqreg_parallel_min_eval(comm = 1L))
      parallel.cond <- FALSE
  }
  if (isTRUE(parallel.cond))
    on.exit(.npRmpi_npqreg_reset_worker_comm_state(comm = 1L), add = TRUE)
  if (!missing(bws) &&
      .npRmpi_npqreg_should_localize(bws) &&
      !isTRUE(getOption("npRmpi.local.regression.mode", FALSE)) &&
      !isTRUE(.npRmpi_autodispatch_in_context()) &&
      !isTRUE(parallel.cond))
    return(.npRmpi_npqreg_eval_local_no_dispatch(
      .npRmpi_eval_without_dispatch(match.call(), parent.frame())
    ))
  if (.npRmpi_autodispatch_active() && !isTRUE(parallel.cond))
    return(.npRmpi_autodispatch_call(match.call(), parent.frame()))

  sc <- sys.call()
  sc.names <- names(sc)

  ## here we check to see if the function was called with tdat =
  ## if it was, we need to catch that and map it to dat =
  ## otherwise the call is passed unadulterated to npudensbw

  bws.named <- any(sc.names == "bws")
  txdat.named <- any(sc.names == "txdat")
  tydat.named <- any(sc.names == "tydat")

  no.bws <- missing(bws)
  no.txdat <- missing(txdat)
  no.tydat <- missing(tydat)
  has.explicit.bws <- (!no.bws) && isa(bws, "condbandwidth")
  bws.formula <- (!no.bws) && inherits(bws, "formula")

  ## if bws was passed in explicitly, do not compute bandwidths
    
  if(txdat.named)
    txdat <- toFrame(txdat)

  if(tydat.named)
    tydat <- toFrame(tydat)

  sc.bw <- sc
  if (bws.formula) {
    sc.bw$`bws` <- NULL
    bws.named <- FALSE
  }

  sc.bw[[1]] <- quote(npcdistbw)
  sc.bw <- .npqreg_strip_fit_controls_from_bw_call(sc.bw)

  if(bws.named){
    sc.bw$bandwidth.compute <- FALSE
  }

  ostxy <- c('txdat','tydat')
  nstxy <- c('xdat','ydat')
  
  m.txy <- match(ostxy, names(sc.bw), nomatch = 0)

  if(any(m.txy > 0)) {
    names(sc.bw)[m.txy] <- nstxy[m.txy > 0]
  }

  use.outer.bandwidth.progress <- !.np_bw_call_uses_nomad_degree_search(
    sc.bw,
    caller_env = parent.frame()
  )

  tbw <- if (!has.explicit.bws) {
    if (use.outer.bandwidth.progress) {
      .np_progress_select_bandwidth_enhanced(
        "Selecting conditional distribution bandwidth",
        .np_eval_bw_call(sc.bw, caller_env = parent.frame())
      )
    } else {
      .np_eval_bw_call(sc.bw, caller_env = parent.frame())
    }
  } else {
    .np_eval_bw_call(sc.bw, caller_env = parent.frame())
  }

  call.args <- list(bws = tbw)
  if (no.bws) {
    call.args$txdat <- txdat
    call.args$tydat <- tydat
  } else {
    if (txdat.named) call.args$txdat <- txdat
    if (tydat.named) call.args$tydat <- tydat
    if ((!bws.named) && (!txdat.named) && (!no.tydat) && (!tydat.named)) {
      call.args <- c(call.args, list(tydat))
    }
  }
  dots <- list(...)
  if (has.explicit.bws)
    fit.dots <- .npqreg_fit_dots(dots)
  else
    fit.dots <- .npqreg_fit_dots(dots, allow.bandwidth.controls = TRUE)
  do.call(npqreg, c(call.args, fit.dots))
}
