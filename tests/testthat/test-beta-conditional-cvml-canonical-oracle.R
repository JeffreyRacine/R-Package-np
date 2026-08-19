beta_cvml_guard <- function(fit) {
  vapply(fit, function(value) {
    if (value > .Machine$double.xmin)
      return(-log(value))
    if (value < -.Machine$double.xmin)
      return(log(-value) - 2 * log(.Machine$double.xmin))
    -log(.Machine$double.xmin)
  }, numeric(1L))
}

beta_cvml_weights <- function(data, bandwidth, type, kernel, order,
                              density = FALSE, evaluation = NULL) {
  arguments <- list(
    bws = bandwidth,
    txdat = data,
    bwtype = type,
    ckertype = kernel,
    ckerorder = order
  )
  if (!is.null(evaluation))
    arguments$exdat <- evaluation
  if (identical(kernel, "beta")) {
    arguments$ckerbound <- "fixed"
    arguments$ckerlb <- 0
    arguments$ckerub <- 1
  }
  if (!density) {
    arguments$return.kernel.weights <- TRUE
    return(do.call(npksum, arguments)$kw)
  }
  # Identity response columns expose the exact response-density operator.
  # For peer kernels public returned `kw` deliberately omits this division.
  arguments$tydat <- diag(nrow(data))
  arguments$bandwidth.divide <- TRUE
  out <- do.call(npksum, arguments)$ksum
  if (is.null(dim(out))) as.numeric(out) else t(out)
}

beta_cvml_adaptive_full_sample_radius <- function(train, k) {
  vapply(seq_along(train), function(index) {
    distance <- abs(train - train[[index]])
    duplicates <- sum(distance == 0) - 1L
    positive <- sort(distance[distance > 0])
    if (duplicates >= k) positive[[1L]] else
      positive[[max(1L, k - duplicates)]]
  }, numeric(1L))
}

beta_cvml_ratio_oracle <- function(x, y, bandwidth, type,
                                    xkernel, ykernel, order) {
  xweights <- beta_cvml_weights(
    x, bandwidth[[2L]], type, xkernel, order, density = FALSE
  )
  if (identical(type, "adaptive_nn") &&
      !identical(xkernel, "beta")) {
    xbandwidth <- beta_cvml_adaptive_full_sample_radius(
      x[[1L]], bandwidth[[2L]]
    )
    xweights <- xweights /
      matrix(xbandwidth, nrow = nrow(x), ncol = nrow(x))
  }
  yweights <- beta_cvml_weights(
    y, bandwidth[[1L]], type, ykernel, order, density = TRUE
  )
  fit <- vapply(seq_len(nrow(x)), function(evaluation) {
    xrow <- xweights[, evaluation]
    yrow <- yweights[, evaluation]
    xrow[[evaluation]] <- 0
    yrow[[evaluation]] <- 0
    sum(xrow * yrow) / sum(xrow)
  }, numeric(1L))

  # The public eval-only wrapper reports the maximized sign of the native
  # negative-log CVML loss.
  -sum(beta_cvml_guard(fit))
}

beta_cvml_native_objective <- function(x, y, bandwidth, type,
                                       xkernel, ykernel, order, compress) {
  old <- options(np.categorical.compress = compress)
  on.exit(options(old), add = TRUE)
  arguments <- list(
    xdat = x,
    ydat = y,
    bws = bandwidth,
    bandwidth.compute = FALSE,
    bwtype = type,
    bwmethod = "cv.ml",
    regtype = "lc",
    cxkertype = xkernel,
    cxkerorder = order,
    cykertype = ykernel,
    cykerorder = order
  )
  if (identical(xkernel, "beta")) {
    arguments$cxkerbound <- "fixed"
    arguments$cxkerlb <- 0
    arguments$cxkerub <- 1
  }
  if (identical(ykernel, "beta")) {
    arguments$cykerbound <- "fixed"
    arguments$cykerlb <- 0
    arguments$cykerub <- 1
  }
  bandwidth_object <- do.call(npcdensbw, arguments)
  evaluator <- getFromNamespace(".npcdensbw_eval_only", "np")
  as.numeric(evaluator(x, y, bandwidth_object)$objective[[1L]])
}

test_that("conditional beta CVML equals the public signed weight-ratio oracle", {
  set.seed(20260801L)
  n <- 31L
  xbeta <- data.frame(x = sort(runif(n, 0.03, 0.97)))
  ybeta <- data.frame(y = sort(runif(n, 0.04, 0.96)))
  xlegacy <- data.frame(x = qnorm(xbeta$x))
  ylegacy <- data.frame(y = qnorm(ybeta$y))

  for (type in c("fixed", "generalized_nn", "adaptive_nn")) {
    bandwidth <- if (identical(type, "fixed")) c(0.19, 0.17) else c(8, 7)
    for (order in c(2L, 4L, 6L, 8L)) {
      for (sides in c("beta_beta", "beta_legacy", "legacy_beta")) {
        xkernel <- if (identical(sides, "legacy_beta")) "gaussian" else "beta"
        ykernel <- if (identical(sides, "beta_legacy")) "gaussian" else "beta"
        x <- if (identical(xkernel, "beta")) xbeta else xlegacy
        y <- if (identical(ykernel, "beta")) ybeta else ylegacy
        expected <- beta_cvml_ratio_oracle(
          x, y, bandwidth, type, xkernel, ykernel, order
        )

        for (compress in c(FALSE, TRUE)) {
          actual <- beta_cvml_native_objective(
            x, y, bandwidth, type, xkernel, ykernel, order, compress
          )
          expect_equal(
            actual,
            expected,
            tolerance = 8e-12,
            info = paste(type, order, sides, "compress", compress)
          )
        }
      }
    }
  }
})
