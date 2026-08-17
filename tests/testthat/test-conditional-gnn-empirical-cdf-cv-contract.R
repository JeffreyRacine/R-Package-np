conditional_gnn_empirical_cdf_radius <- function(train, evaluation, k,
                                                  excluded) {
  keep <- setdiff(seq_along(train), excluded)
  sort(abs(train[keep] - evaluation), method = "radix")[[k]]
}

conditional_gnn_empirical_cdf_oracle <- function(x, y, kx, ky) {
  x <- as.matrix(x)
  y <- as.matrix(y)
  n <- nrow(x)
  total <- 0

  for (i in seq_len(n)) {
    keep <- setdiff(seq_len(n), i)
    hx <- vapply(seq_len(ncol(x)), function(d) {
      conditional_gnn_empirical_cdf_radius(
        x[, d], x[i, d], kx[[d]], i
      )
    }, numeric(1L))
    xrow <- rep.int(1, length(keep))
    for (d in seq_len(ncol(x)))
      xrow <- xrow * dnorm((x[i, d] - x[keep, d]) / hx[[d]])
    xrow <- xrow / sum(xrow)

    for (j in seq_len(n)) {
      if (i == j)
        next
      hy <- vapply(seq_len(ncol(y)), function(d) {
        conditional_gnn_empirical_cdf_radius(
          y[, d], y[j, d], ky[[d]], c(i, j)
        )
      }, numeric(1L))
      yrow <- rep.int(1, length(keep))
      for (d in seq_len(ncol(y)))
        yrow <- yrow * pnorm((y[j, d] - y[keep, d]) / hy[[d]])
      fit <- sum(xrow * yrow)
      indicator <- as.integer(all(y[i, ] <= y[j, ]))
      total <- total + (indicator - fit)^2
    }
  }
  total / (n * (n - 1L))
}

conditional_gnn_empirical_cdf_native <- function(x, y, kx, ky, tree) {
  options(np.tree = tree)
  bw <- npcdistbw(
    xdat = as.data.frame(x), ydat = as.data.frame(y),
    bws = c(ky, kx), bandwidth.compute = FALSE,
    bwmethod = "cv.ls", bwtype = "generalized_nn", bwscaling = FALSE,
    cxkertype = "gaussian", cykertype = "gaussian", regtype = "lc"
  )
  np:::.npcdistbw_eval_only(
    as.data.frame(x), as.data.frame(y), bws = bw,
    do.full.integral = TRUE
  )$objective[[1L]]
}

test_that("empirical conditional-CDF CV uses pair-local generalized-NN radii", {
  old <- options(np.messages = FALSE, np.tree = FALSE,
                 np.macMseries.accelerate = FALSE)
  on.exit(options(old), add = TRUE)
  set.seed(20260822)
  n <- 9L
  x <- data.frame(x = sort(runif(n, -0.9, 1.1)))
  y <- data.frame(y = 0.55 * x$x + rnorm(n, sd = 0.23))
  expected <- conditional_gnn_empirical_cdf_oracle(x, y, 3L, 4L)
  observed <- vapply(c(FALSE, TRUE), function(tree) {
    conditional_gnn_empirical_cdf_native(x, y, 3L, 4L, tree)
  }, numeric(1L))

  expect_equal(observed, rep(expected, 2L), tolerance = 5e-10)
  expect_identical(observed[[1L]], observed[[2L]])
})

test_that("empirical conditional-CDF CV composes response-coordinate exclusions", {
  old <- options(np.messages = FALSE, np.tree = FALSE,
                 np.macMseries.accelerate = FALSE)
  on.exit(options(old), add = TRUE)
  set.seed(20260823)
  n <- 9L
  x <- data.frame(x = sort(runif(n, -0.7, 0.8)))
  y <- data.frame(
    y1 = 0.4 * x$x + rnorm(n, sd = 0.18),
    y2 = -0.2 * x$x + rnorm(n, sd = 0.21)
  )
  expected <- conditional_gnn_empirical_cdf_oracle(x, y, 3L, c(3L, 4L))
  observed <- vapply(c(FALSE, TRUE), function(tree) {
    conditional_gnn_empirical_cdf_native(x, y, 3L, c(3L, 4L), tree)
  }, numeric(1L))

  expect_equal(observed, rep(expected, 2L), tolerance = 8e-10)
  expect_identical(observed[[1L]], observed[[2L]])
})

test_that("empirical conditional-CDF CV preserves occurrences and the last admissible k", {
  old <- options(np.messages = FALSE, np.tree = FALSE,
                 np.macMseries.accelerate = FALSE)
  on.exit(options(old), add = TRUE)

  x <- data.frame(x = c(-0.8, -0.8, -0.25, 0, 0,
                        0.31, 0.55, 0.55, 0.82, 1.05))
  y <- data.frame(y = c(-0.4, -0.4, -0.08, 0.12, 0.12,
                        0.34, 0.48, 0.48, 0.73, 0.9))
  expected.duplicates <- conditional_gnn_empirical_cdf_oracle(
    x, y, 3L, 3L
  )
  observed.duplicates <- vapply(c(FALSE, TRUE), function(tree) {
    conditional_gnn_empirical_cdf_native(x, y, 3L, 3L, tree)
  }, numeric(1L))

  set.seed(20260825)
  n <- 10L
  x <- data.frame(x = sort(runif(n, -0.8, 1)))
  y <- data.frame(y = 0.6 * x$x + rnorm(n, sd = 0.24))
  saturated.k <- n - 2L
  expected.saturated <- conditional_gnn_empirical_cdf_oracle(
    x, y, 3L, saturated.k
  )
  observed.saturated <- vapply(c(FALSE, TRUE), function(tree) {
    conditional_gnn_empirical_cdf_native(x, y, 3L, saturated.k, tree)
  }, numeric(1L))

  expect_equal(observed.duplicates, rep(expected.duplicates, 2L),
               tolerance = 2e-10)
  expect_identical(observed.duplicates[[1L]], observed.duplicates[[2L]])
  expect_equal(observed.saturated, rep(expected.saturated, 2L),
               tolerance = 2e-9)
  expect_identical(observed.saturated[[1L]], observed.saturated[[2L]])
})

test_that("empirical conditional-CDF CV counts only off-diagonal pairs", {
  old <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old), add = TRUE)
  set.seed(20260824)
  n <- 8L
  x <- data.frame(x = runif(n))
  y <- data.frame(y = 0.3 * x$x + rnorm(n, sd = 0.2))
  bw <- npcdistbw(
    xdat = x, ydat = y, bws = c(0.31, 0.43),
    bandwidth.compute = FALSE, bwmethod = "cv.ls",
    bwtype = "fixed", bwscaling = FALSE,
    cxkertype = "gaussian", cykertype = "gaussian", regtype = "lc"
  )
  objective <- np:::.npcdistbw_eval_only(
    x, y, bws = bw, do.full.integral = TRUE
  )$objective[[1L]]

  total <- 0
  for (i in seq_len(n)) {
    keep <- setdiff(seq_len(n), i)
    xrow <- dnorm((x$x[[i]] - x$x[keep]) / 0.43)
    xrow <- xrow / sum(xrow)
    for (j in seq_len(n)) {
      if (i == j)
        next
      fit <- sum(xrow * pnorm((y$y[[j]] - y$y[keep]) / 0.31))
      total <- total + (as.integer(y$y[[i]] <= y$y[[j]]) - fit)^2
    }
  }
  expected <- total / (n * (n - 1L))
  expect_equal(objective, expected, tolerance = 5e-10)
})
