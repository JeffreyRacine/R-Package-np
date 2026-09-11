# One fixed-coordinate density-ratio uncertainty owner, shared by fitting,
# prediction and asymptotic plot payloads. No work is done for se = FALSE.
.npcopula_density_se_local <- function(bws, data, xgrid, progress.context = NULL) {
  if (!identical(bws[["type", exact = TRUE]], "fixed"))
    stop("internal copula influence requires fixed bandwidths", call. = FALSE)
  n <- nrow(data)
  m <- nrow(xgrid)
  p <- length(bws[["xnames", exact = TRUE]])
  if (n < 1L || p < 1L)
    stop("internal copula influence requires nonempty training data", call. = FALSE)
  if (m == 0L)
    return(numeric(0L))

  # Account for all p resident marginal tiles, not just a single matrix.
  tile <- max(1L, min(64L, floor(1048576 / (as.double(n) * p))))
  marginal.bws <- lapply(seq_len(p), function(j)
    .npcopula_marginal_bw(bws, data, j, target = "density", kbandwidth = TRUE))
  training <- lapply(seq_len(p), function(j)
    .npcopula_marginal_data(bws, data, j))
  evaluation <- lapply(seq_len(p), function(j)
    .npcopula_marginal_eval_data(bws, xgrid, j))
  out <- numeric(m)
  for (start in seq.int(1L, m, by = tile)) {
    rows <- seq.int(start, min(m, start + tile - 1L))
    weights <- lapply(seq_len(p), function(j) {
      npksum(bws = marginal.bws[[j]], txdat = training[[j]],
             exdat = evaluation[[j]][rows, , drop = FALSE],
             operator = "normal", bandwidth.divide = FALSE,
             return.kernel.weights = TRUE)[["kw", exact = TRUE]]
    })
    out[rows] <- .Call("C_np_copula_density_se", weights, PACKAGE = "np")
    weights <- NULL
    if (!is.null(progress.context))
      progress.context$state <- .np_progress_step_at(
        progress.context$state, now = .np_progress_now(),
        done = if (isTRUE(progress.context$use.bootstrap.done)) max(rows) else
          progress.context$done,
        detail = sprintf("asymptotic standard errors (%d/%d rows)", max(rows), m))
  }
  out
}

.npcopula_density_se_fixed <- function(bws, data, xgrid, progress.context = NULL) {
  if (is.null(progress.context)) {
    progress.context <- new.env(parent = emptyenv())
    progress.context$state <- .np_progress_begin(
      "Copula density standard errors", total = nrow(xgrid), surface = "copula")
    progress.context$use.bootstrap.done <- TRUE
    on.exit(.np_progress_end(progress.context$state), add = TRUE)
  }
  out <- .npcopula_density_se_local(bws, data, xgrid, progress.context)
  if (anyNA(out))
    .np_warning(paste0(
      "Copula density standard errors are unavailable for ",
      sum(is.na(out)), " evaluation row(s): a marginal density is zero ",
      "or the ratio influence is nonfinite. These standard errors are NA."
    ), call. = FALSE)
  out
}
