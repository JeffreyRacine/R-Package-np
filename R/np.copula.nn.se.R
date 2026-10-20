# Pointwise NN density-ratio influence at held physical query coordinates.
# The unchanged density primitives own kernels, ANN radii and support faces.
# Factors returned by those primitives are already centered and scaled.
.npcopula_density_se_nn_norm <- function(score) {
  if (any(!is.finite(score))) return(NA_real_)
  scale <- max(abs(score))
  if (scale == 0) return(0)
  z <- score / scale
  z <- z - mean(z)
  value <- scale * sqrt(sum(z * z))
  if (is.finite(value)) value else NA_real_
}

.npcopula_density_se_nn_factor <- function(bws, data, evaluation) {
  tryCatch(.np_plot_density_unconditional(bws, data, evaluation),
    np_plot_joint_unavailable = function(e) {
      if (conditionMessage(e) %in% c("nonregular ANN spacing geometry",
                                     "nonregular ANN support face")) return(NULL)
      stop(e)
    })
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
    .npcopula_marginal_bw(bws, data, j, target = "density"))
  training <- lapply(seq_len(p), function(j)
    .npcopula_marginal_data(bws, data, j))
  evaluation <- lapply(seq_len(p), function(j)
    .npcopula_marginal_eval_data(bws, xgrid, j))
  if (p == 1L &&
      (!identical(data[[1L]], training[[1L]][[1L]]) ||
       !identical(xgrid[[1L]], evaluation[[1L]][[1L]]) ||
       !identical(.np_make_kbandwidth_unconditional(bws, data),
                  .np_make_kbandwidth_unconditional(marginal.bws[[1L]], training[[1L]]))))
    stop("NN copula marginal primitive identity failed", call. = FALSE)
  # Conservative allowance for joint/child factors, derivatives, faces and R
  # copies. Geometry is prepared by the existing primitive per tile invocation.
  tile <- max(1L, min(64L, floor(4194304 / (8 * as.double(n) * (80 + 12 * p)))))
  # The existing face primitive reports a whole-call regularity condition.
  # Query-local calls preserve row-local availability without a retry/fallback.
  if (ann && identical(bws[["ckertype", exact = TRUE]], "uniform")) tile <- 1L
  out <- rep.int(NA_real_, m)
  for (start in seq.int(1L, m, by = tile)) {
    rows <- seq.int(start, min(m, start + tile - 1L))
    joint <- .npcopula_density_se_nn_factor(bws, data, xgrid[rows, , drop = FALSE])
    if (is.null(joint)) next
    if (p == 1L) {
      # Reuse the same primitive; no subtraction of approximately equal factors.
      valid <- is.finite(joint$estimate) & joint$estimate != 0 &
        apply(is.finite(joint$A), 1L, all)
      out[rows[valid]] <- 0
    } else {
      denominator <- rep.int(1, length(rows))
      marginal <- matrix(0, length(rows), n)
      valid <- rep.int(TRUE, length(rows))
      for (j in seq_len(p)) {
        child <- .npcopula_density_se_nn_factor(marginal.bws[[j]], training[[j]],
          evaluation[[j]][rows, , drop = FALSE])
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
        # No division by J: a zero numerator with defined margins is valid.
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
  out <- .npcopula_density_se_nn_local(bws, data, xgrid, progress.context)
  if (anyNA(out))
    .np_warning(paste0("Copula density standard errors are unavailable for ",
      sum(is.na(out)), " evaluation row(s): a marginal density is zero, ",
      "NN spacing or support-face geometry is nonregular, or the ratio ",
      "influence is nonfinite. These standard errors are NA."), call. = FALSE)
  out
}
