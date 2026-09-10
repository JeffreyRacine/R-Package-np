# Single-index inference uses the regression engine on the fitted scalar index.
# These helpers do not select bandwidths or resample beta.
# Sum response moments and their common denominator without kernel matrices.
# Stabilization belongs to the caller, not to this raw-moment owner.
.np_index_kernel_moments <- function(y, ...) {
  y <- as.matrix(y)
  sums <- .np_index_kernel_sum(
    ..., tydat = rep.int(1.0, nrow(y)), weights = cbind(y, 1.0)
  )$ksum
  list(numerator = sums[seq_len(ncol(y)), , drop = FALSE],
       denominator = as.numeric(sums[ncol(y) + 1L, ]))
}

.np_index_kernel_sum <- function(..., bwtype,
                                 bandwidth.divide = identical(bwtype, "adaptive_nn")) {
  if (is.null(.np_progress_runtime$fit_forward))
    return(npksum(..., bwtype = bwtype,
                  bandwidth.divide = bandwidth.divide))
  # The parent reports activity, not a percentage across heterogeneous calls.
  # Enable the existing native row callbacks only during this scoped activity.
  args <- list(...)
  total <- max(NROW(args[["txdat"]]), NROW(args[["exdat"]]))
  .np_with_compiled_fit_progress(
    label = "Fitting single-index model", total = total,
    expr = npksum(..., bwtype = bwtype,
                  bandwidth.divide = bandwidth.divide))
}

# Preserve the raw kernel-sum owner and its arithmetic; only exceptional
# external rows need computed-weight evidence. Training ratios stay strict.
.np_index_normalized_mean <- function(sums, txdat, exdat, bws,
                                       allow.empty.rows = FALSE) {
  denominator <- sums[2L, 2L, ]
  denominator <- .np_normalization_denominator(
    denominator, "npindex", allow.empty.rows = allow.empty.rows,
    zero.rows = .np_indexhat_zero_moment_rows(list(
      txdat = toFrame(txdat), exdat = toFrame(exdat), bws = bws$bw,
      bwtype = bws$type, ckertype = bws$ckertype,
      ckerorder = bws$ckerorder, ckerbound = bws$ckerbound), denominator))
  .np_normalization_finish(sums[1L, 2L, ] / denominator, denominator,
                           "npindex", defer.empty.rows = TRUE)
}

# The LC regression owner predates explicit empty-row metadata. Retain that
# owner and collect evidence only when it actually returned missing means.
.np_index_fit_rows <- function(fit, txdat, exdat, bws, allow.empty.rows) {
  if (!isTRUE(allow.empty.rows) ||
      !is.null(attr(fit, ".np.empty.rows", exact = TRUE)) || !anyNA(fit$mean))
    return(fit)
  rows <- which(is.na(fit$mean))
  args <- list(txdat = toFrame(txdat),
    exdat = toFrame(exdat)[rows, , drop = FALSE], bws = bws$bw,
    bwtype = bws$type, ckertype = bws$ckertype,
    ckerorder = bws$ckerorder, ckerbound = bws$ckerbound)
  moments <- do.call(.np_index_kernel_moments,
    c(list(y = rep.int(1.0, NROW(txdat))), args))
  empty <- .np_indexhat_zero_moment_rows(args, moments$denominator)
  if (any(empty)) {
    flags <- integer(length(fit$mean))
    flags[rows[empty]] <- 1L
    attr(fit, ".np.empty.rows") <- flags
  }
  fit
}

.np_index_asymptotic_outputs <- function(fit, beta, gradients = FALSE) {
  out <- list(merr = as.double(fit$merr))
  if (gradients) {
    out$gerr <- as.vector(fit$gerr[, 1L]) %o% abs(as.vector(beta))
  }
  out
}

.np_index_covariance_inverse <- function(information) {
  factor <- tryCatch(chol(information), error = function(e) {
    stop(paste0("npindex(): asymptotic coefficient covariance could not be computed: ",
                "the free-coefficient information matrix is singular or not positive definite. ",
                "Use se = FALSE for point estimates without inference. ",
                "Original factorization error: ", conditionMessage(e)),
         call. = FALSE)
  })
  chol2inv(factor)
}

.np_index_refit_hint <- function(expr, gradients = FALSE, se = FALSE) {
  .np_se_refit_hint(
    expr, "npindex",
    switches = paste0("gradients = ", if (gradients) "TRUE" else "FALSE",
                      ", se = ", if (se) "TRUE" else "FALSE"))
}
