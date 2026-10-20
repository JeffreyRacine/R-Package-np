# Rank-local pointwise NN density-ratio uncertainty. Existing package-local
# primitives retain ownership of kernels, donor radii and uniform support faces.
.npcopula_density_se_nn_norm <- function(score) {
  if (any(!is.finite(score))) return(NA_real_)
  scale <- max(abs(score))
  if (scale == 0) return(0)
  z <- score / scale
  z <- z - mean(z)
  value <- scale * sqrt(sum(z * z))
  if (is.finite(value)) value else NA_real_
}

# All inputs here are retained bandwidth metadata, not a fresh kernel-order
# request. Reuse the common exact advisory policy; other conditions propagate.
.npcopula_density_se_nn_retained <- function(bws, expr) {
  .np_retained_uniform_constructor(function(ckertype) force(expr),
    bws[["ckertype", exact = TRUE]])
}

.npcopula_density_se_nn_factor <- function(bws, data, evaluation, prepared = NULL) {
  if (isTRUE(prepared[["unavailable", exact = TRUE]])) return(NULL)
  .npcopula_density_se_nn_retained(bws,
    tryCatch(if (is.null(prepared))
      .np_plot_density_unconditional(bws, data, evaluation) else
      .np_plot_density_ann(bws, data, evaluation, TRUE, prepared),
    np_plot_joint_unavailable = function(e) {
      if (conditionMessage(e) %in% c("nonregular ANN spacing geometry",
                                     "nonregular ANN support face")) return(NULL)
      stop(e)
    }))
}

.npcopula_density_se_nn_local <- function(bws, data, xgrid,
                                         progress.context = NULL) {
  n <- nrow(data); m <- nrow(xgrid); p <- ncol(data)
  if (!bws[["type", exact = TRUE]] %in% c("generalized_nn", "adaptive_nn") ||
      n < 1L || p < 1L)
    stop("invalid NN copula uncertainty input", call. = FALSE)
  if (!m) return(numeric(0L))
  attr(data, "na.action") <- NULL
  attr(xgrid, "na.action") <- NULL
  ann <- identical(bws[["type", exact = TRUE]], "adaptive_nn") &&
    any(bws[["icon", exact = TRUE]])
  marginal.bws <- lapply(seq_len(p), function(j)
    .npcopula_density_se_nn_retained(bws,
      .npcopula_marginal_bw(bws, data, j, target = "density")))
  training <- lapply(seq_len(p), function(j)
    .npcopula_marginal_data(bws, data, j))
  evaluation <- lapply(seq_len(p), function(j)
    .npcopula_marginal_eval_data(bws, xgrid, j))
  if (p == 1L) {
    joint.spec <- .npcopula_density_se_nn_retained(bws,
      .np_make_kbandwidth_unconditional(bws, data))
    marginal.spec <- .npcopula_density_se_nn_retained(marginal.bws[[1L]],
      .np_make_kbandwidth_unconditional(marginal.bws[[1L]], training[[1L]]))
    if (joint.spec[["ncon", exact = TRUE]] == 0L &&
        marginal.spec[["ncon", exact = TRUE]] == 0L) {
      # Continuous metadata cannot distinguish two purely categorical kernels.
      inactive <- c("ckertype", "ckerorder", "ckerbound", "ckerlb", "ckerub", "pckertype")
      joint.spec[inactive] <- NULL
      marginal.spec[inactive] <- NULL
    }
    if (!identical(data[[1L]], training[[1L]][[1L]]) ||
        !identical(xgrid[[1L]], evaluation[[1L]][[1L]]) ||
        !identical(joint.spec, marginal.spec))
      stop("NN copula marginal primitive identity failed", call. = FALSE)
  }
  # Contexts belong only to this local invocation. Reuse training geometry,
  # never query-specific kernels, support faces or influence factors.
  joint.prepared <- NULL
  marginal.prepared <- vector("list", p)
  if (ann && identical(bws[["ckertype", exact = TRUE]], "uniform")) {
    prepare <- function(b, x) {
      if (!any(b[["icon", exact = TRUE]])) return(NULL)
      .npcopula_density_se_nn_retained(b,
        tryCatch(.np_plot_density_ann_prepare(b, x, TRUE, uniform.faces = TRUE),
          np_plot_joint_unavailable = function(e) {
            if (identical(conditionMessage(e), "nonregular ANN spacing geometry"))
              return(list(unavailable = TRUE))
            stop(e)
          }))
    }
    joint.prepared <- prepare(bws, data)
    if (p > 1L && !isTRUE(joint.prepared[["unavailable", exact = TRUE]]))
      marginal.prepared <- lapply(seq_len(p), function(j)
        prepare(marginal.bws[[j]], training[[j]]))
  }
  tile <- max(1L, min(64L, floor(4194304 / (8 * as.double(n) * (80 + 12 * p)))))
  # Face availability is query-local; the existing primitive signals a whole
  # input condition. Single-query calls preserve every other regular row.
  if (ann && identical(bws[["ckertype", exact = TRUE]], "uniform")) tile <- 1L
  out <- rep.int(NA_real_, m)
  for (start in seq.int(1L, m, by = tile)) {
    rows <- seq.int(start, min(m, start + tile - 1L))
    joint <- .npcopula_density_se_nn_factor(bws, data, xgrid[rows, , drop = FALSE],
      joint.prepared)
    if (is.null(joint)) next
    if (p == 1L) {
      valid <- is.finite(joint$estimate) & joint$estimate != 0 &
        apply(is.finite(joint$A), 1L, all)
      out[rows[valid]] <- 0
    } else {
      denominator <- rep.int(1, length(rows))
      marginal <- matrix(0, length(rows), n)
      valid <- rep.int(TRUE, length(rows))
      for (j in seq_len(p)) {
        child <- .npcopula_density_se_nn_factor(marginal.bws[[j]], training[[j]],
          evaluation[[j]][rows, , drop = FALSE], marginal.prepared[[j]])
        if (is.null(child)) { valid[] <- FALSE; break }
        if (ann && !any(marginal.bws[[j]][["icon", exact = TRUE]]))
          child$A <- child$A * sqrt(n / (n - 1))
        valid <- valid & is.finite(child$estimate) & child$estimate != 0
        denominator <- denominator * child$estimate
        marginal <- marginal + child$A / child$estimate
      }
      valid <- valid & is.finite(denominator) & denominator != 0 &
        is.finite(joint$estimate)
      for (q in which(valid)) {
        score <- joint$A[q, ] / denominator[q] -
          (joint$estimate[q] / denominator[q]) * marginal[q, ]
        out[rows[q]] <- .npcopula_density_se_nn_norm(score)
      }
    }
    if (!is.null(progress.context))
      progress.context$state <- .np_progress_step_at(
        progress.context$state, now = .np_progress_now(),
        done = if (isTRUE(progress.context$use.bootstrap.done)) max(rows) else
          progress.context$done,
        detail = sprintf("asymptotic standard errors (%d/%d rows)", max(rows), m))
  }
  out
}

.npcopula_density_se_nn <- function(bws, data, xgrid, progress.context = NULL) {
  ann <- identical(bws[["type", exact = TRUE]], "adaptive_nn") &&
    any(bws[["icon", exact = TRUE]])
  if (ann) {
    reason <- .np_ann_se_domain(bws, nrow(data), TRUE)
    if (!is.null(reason)) {
      .np_warning(paste0("Adaptive-NN standard errors are unavailable: ", reason,
        ". Point estimates are retained; standard errors are NA."), call. = FALSE)
      return(rep.int(NA_real_, nrow(xgrid)))
    }
  }
  if (is.null(progress.context)) {
    progress.context <- new.env(parent = emptyenv())
    progress.context$state <- .np_progress_begin(
      "Copula density standard errors", total = nrow(xgrid), surface = "copula")
    progress.context$use.bootstrap.done <- TRUE
    on.exit(.np_progress_end(progress.context$state), add = TRUE)
  }
  out <- .npRmpi_copula_density_se_nn(bws, data, xgrid, progress.context)
  if (anyNA(out))
    .np_warning(paste0("Copula density standard errors are unavailable for ",
      sum(is.na(out)), " evaluation row(s): a marginal density is zero, ",
      "NN spacing or support-face geometry is nonregular, or the ratio ",
      "influence is nonfinite. These standard errors are NA."), call. = FALSE)
  out
}

# The established fixed-copula transport is the analogue, kept unchanged.
# Complete rows stay rank-local; no influence matrix is gathered or leased.
.npRmpi_copula_density_se_nn <- function(bws, data, xgrid, progress.context = NULL) {
  .npRmpi_require_active_slave_pool(where = "npcopula standard errors")
  if (isTRUE(getOption("npRmpi.local.regression.mode", FALSE)))
    return(.npcopula_density_se_nn_local(bws, data, xgrid, progress.context))
  m <- nrow(xgrid)
  if (!m) return(numeric(0L))
  if (.npRmpi_autodispatch_called_from_bcast()) {
    rank <- as.integer(mpi.comm.rank(comm = 1L))
    size <- as.integer(mpi.comm.size(comm = 1L))
    base <- m %/% size; extra <- m %% size
    count <- base + as.integer(rank < extra)
    start <- rank * base + min(rank, extra) + 1L
    rows <- if (count > 0L) seq.int(start, length.out = count) else integer(0L)
    local <- .npRmpi_capture_local_work(.npRmpi_with_local_regression(
      .npcopula_density_se_nn_local(bws, data, xgrid[rows, , drop = FALSE], progress.context)))
    parts <- mpi.allgather.Robj(list(list(rows = rows, value = local)), comm = 1L)
    out <- numeric(m); seen <- integer(m)
    for (wrapped in parts) {
      part <- if (is.list(wrapped) && length(wrapped) == 1L &&
                  is.list(wrapped[[1L]])) wrapped[[1L]] else wrapped
      if (!is.list(part) || is.null(part[["rows"]]) || is.null(part[["value"]]))
        stop("NN copula SE gathered a malformed result", call. = FALSE)
      if (inherits(part[["value"]], "npRmpi_local_failure"))
        .npRmpi_raise_completed_failure(part[["value"]])
      ids <- part[["rows"]]
      if (!is.integer(ids) || anyNA(ids) || any(ids < 1L | ids > m) ||
          anyDuplicated(ids) || !is.double(part[["value"]]) ||
          length(part[["value"]]) != length(ids))
        stop("NN copula SE gathered invalid row ownership", call. = FALSE)
      out[ids] <- part[["value"]]; seen[ids] <- seen[ids] + 1L
    }
    if (any(seen != 1L))
      stop("NN copula SE row ownership was incomplete or duplicated", call. = FALSE)
    return(out)
  }
  tasks <- .npRmpi_bootstrap_chunk_tasks(B = m, chunk.size = 64L, with.seeds = FALSE)
  snapshot <- .npRmpi_autodispatch_option_snapshot()
  worker <- function(task, bws, data, xgrid, snapshot) {
    old <- options(snapshot)
    on.exit(options(old), add = TRUE)
    rows <- seq.int(task$start, length.out = task$bsz)
    values <- .npRmpi_with_local_regression(
      .npcopula_density_se_nn_local(bws, data, xgrid[rows, , drop = FALSE]))
    matrix(values, ncol = 1L)
  }
  seed.state <- list(
    exists.seed = exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE),
    save.seed = get0(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
  on.exit(.np_seed_exit(seed.state, remove_if_absent = TRUE), add = TRUE)
  result <- .npRmpi_bootstrap_run_fanout(
    progress.unit = "row", tasks = tasks, worker = worker, ncol.out = 1L,
    what = "npcopula-density-nn-se", progress.label = "Copula density standard errors",
    profile.where = "npcopula-density-nn-se", comm = 1L,
    prefer.local.single_worker = FALSE, master_local_chunk = TRUE,
    bws = bws, data = data, xgrid = xgrid, snapshot = snapshot,
    progress.context = progress.context)
  as.numeric(result[, 1L])
}
