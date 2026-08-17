conditional_gnn_radius <- function(train, evaluation, k, excluded) {
  keep <- setdiff(seq_along(train), excluded)
  sort(abs(train[keep] - evaluation), method = "radix")[[k]]
}

conditional_gnn_cvml_oracle <- function(x, y, kx, ky, degree = 0L) {
  x <- as.matrix(x)
  y <- as.numeric(y)
  n <- nrow(x)

  fitted <- vapply(seq_len(n), function(i) {
    hx <- vapply(seq_len(ncol(x)), function(l) {
      conditional_gnn_radius(x[, l], x[i, l], kx[[l]], i)
    }, numeric(1L))
    hy <- conditional_gnn_radius(y, y[[i]], ky, i)
    wx <- rep.int(1, n)
    for (l in seq_len(ncol(x)))
      wx <- wx * dnorm((x[i, l] - x[, l]) / hx[[l]]) / hx[[l]]
    wy <- dnorm((y[[i]] - y) / hy) / hy
    wx[[i]] <- 0
    wy[[i]] <- 0

    if (degree == 0L) {
      influence <- wx / sum(wx)
    } else {
      design <- cbind(1, x[, 1L] - x[i, 1L])
      gram <- crossprod(design, wx * design)
      influence <- drop(c(1, 0) %*%
        solve(gram, t(design) * rep(wx, each = 2L)))
    }
    sum(influence * wy)
  }, numeric(1L))

  stopifnot(all(is.finite(fitted)), all(fitted > 0))
  sum(log(fitted))
}

conditional_gnn_cvml_bw <- function(x, y, kx, ky, degree) {
  arguments <- list(
    xdat = as.data.frame(x), ydat = data.frame(y = y),
    bws = c(ky, kx), bandwidth.compute = FALSE,
    bwmethod = "cv.ml", bwtype = "generalized_nn", bwscaling = FALSE,
    cxkertype = "gaussian", cykertype = "gaussian",
    regtype = if (degree == 0L) "lc" else "lp"
  )
  if (degree > 0L) {
    arguments$basis <- "glp"
    arguments$degree <- rep.int(degree, ncol(as.matrix(x)))
  }
  do.call(npcdensbw, arguments)
}

test_that("conditional generalized-NN CVML deletes focal radius occurrences", {
  old <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old), add = TRUE)
  set.seed(20260816)
  n <- 13L
  x1 <- sort(runif(n, -0.9, 1.1))
  x2 <- runif(n, -0.6, 0.8)
  y <- 0.8 * x1 - 0.35 * x2 + rnorm(n, sd = 0.17)
  cases <- list(
    lc1 = list(x = data.frame(x1 = x1), kx = 4L, ky = 5L, degree = 0L),
    lc2 = list(x = data.frame(x1 = x1, x2 = x2), kx = c(4L, 6L), ky = 5L, degree = 0L),
    ll1 = list(x = data.frame(x1 = x1), kx = 6L, ky = 5L, degree = 1L)
  )

  for (case in cases) {
    expected <- conditional_gnn_cvml_oracle(
      case$x, y, case$kx, case$ky, case$degree
    )
    bw <- conditional_gnn_cvml_bw(
      case$x, y, case$kx, case$ky, case$degree
    )
    for (tree in c(FALSE, TRUE)) {
      options(np.tree = tree)
      observed <- np:::.npcdensbw_eval_only(
        case$x, data.frame(y = y), bw
      )$objective[[1L]]
      expect_equal(observed, expected, tolerance = 2e-10)
    }
  }
})

test_that("conditional generalized-NN CVML rejects zero literal radii", {
  old <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old), add = TRUE)
  x <- data.frame(x = c(0, 0, 0, 1, 2))
  y <- data.frame(y = c(0, 0, 0, 0.5, 1))
  bw <- npcdensbw(
    xdat = x, ydat = y, bws = c(2L, 2L), bandwidth.compute = FALSE,
    bwmethod = "cv.ml", bwtype = "generalized_nn", bwscaling = FALSE
  )

  observed <- np:::.npcdensbw_eval_only(
    x, y, bw, invalid.penalty = "dbmax"
  )$objective[[1L]]
  expect_identical(observed, -.Machine$double.xmax)
})
