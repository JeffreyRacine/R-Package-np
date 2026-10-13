# Bounded adapter contracts; full release/suite work belongs to the release lane.
test_that("four-family calibration is scale-safe and refuses false degeneracy", {
  A <- rbind(c(1, -2, 1, 0), c(2, 0, -1, -1), rep(0, 4))
  calibrate <- function(a, s = c(1, 2, 0)) .np_with_seed(41,
    .np_plot_four_calibrate(a, s, .05, 9999L, "test"))
  reference <- calibrate(A)
  expect_identical(reference$critical, calibrate(A * c(1e-250, 1e250, 1))$critical)
  expect_error(calibrate(A * 1e-250, c(0, 2, 0)), "disagree about zero")
  expect_error(calibrate(matrix(1e308, 2, 4), c(1, 1)), "not representable")
  expect_equal(reference$grid.size, 3L)
  expect_error(.np_plot_four_budget(1e6, 1e6, 1), "memory budget")
})

test_that("quantile influence excludes support clamping even near the boundary", {
  set.seed(71)
  x <- data.frame(x = seq(-1, 1, length.out = 50))
  y <- .3 * x$x + rnorm(50, sd = .6)
  ex <- data.frame(x = 0)
  b <- npcdistbw(xdat = x, ydat = y, bws = c(.4, .5), bandwidth.compute = FALSE)
  low <- .np_plot_density_conditional(b, x, data.frame(y = y), ex,
    data.frame(y = min(y)), cdf = TRUE)$estimate
  tau <- low - 1e-10
  q <- .npqreg_invert_selected_cdf(bws = b, xdat = x, ydat = data.frame(y = y),
    exdat = ex, tau = tau, tol = 1e-8, small = 1e-5, itmax = 1000)
  expect_equal(as.numeric(q), min(y))
  for (root in list(q, as.numeric(q)))
    expect_error(.np_plot_qreg_influence(b, x, data.frame(y = y), ex, tau,
      root, 1), "interior, unclamped")
  fit <- npqreg(bws = b, txdat = x, tydat = y, exdat = ex,
    tau = low + .01, se = TRUE)
  expect_true(all(is.finite(.np_plot_qreg_influence(b, x, data.frame(y = y),
    ex, low + .01, as.numeric(fitted(fit)), as.numeric(se(fit)))$A)))
})

test_that("quantile raw covariance uses the shared inverse-CDF observation rows", {
  set.seed(31)
  n <- 48L; x <- data.frame(x = seq(-1, 1, length.out = n))
  y <- sin(x$x) + rnorm(n, sd = .3); ex <- x[c(7, 17, 29, 41), , drop = FALSE]
  b <- npcdistbw(xdat = x, ydat = y, bws = c(.25, .4),
    bandwidth.compute = FALSE, regtype = "lc")
  for (degree in 1:2) {
    positive <- npcdistbw(xdat = x, ydat = y, bws = c(.25, .4),
      bandwidth.compute = FALSE, regtype = "lp", degree = degree)
    expect_error(.np_plot_four_conditional_admit(positive, "quantile"),
      "positive-degree owner covariance is not qualified")
  }
  fit <- npqreg(bws = b, txdat = x, tydat = y, exdat = ex, tau = .4, se = TRUE)
  q <- fitted(fit); actual <- .np_plot_qreg_influence(b, x, data.frame(y = y),
    ex, .4, q, as.numeric(se(fit)))
  # Independent monomial least squares; no package basis/influence owner.
  W <- outer(x$x, 0, "^"); E <- outer(ex$x, 0, "^")
  oracle <- matrix(0, nrow(ex), n)
  for (j in seq_len(nrow(ex))) {
    w <- dnorm((x$x - ex$x[j]) / .4)
    M <- crossprod(W, w * W)
    response <- pnorm((q[j] - y) / .25)
    beta <- solve(M, crossprod(W, w * response))
    slope <- sum(E[j, ] * solve(M, crossprod(W, w * dnorm((q[j] - y) / .25) / .25)))
    oracle[j, ] <- -w * drop(W %*% solve(M, E[j, ])) *
      (response - drop(W %*% beta)) / slope
  }
  oracle <- oracle - rowMeans(oracle)
  expect_lt(max(abs(tcrossprod(actual$A) - tcrossprod(oracle))), 2e-11)
  perm <- c(4, 2, 2, 1)
  reordered <- .np_plot_qreg_influence(b, x, data.frame(y = y),
    ex[perm, , drop = FALSE], .4, q[perm], as.numeric(se(fit))[perm])
  expect_lt(max(abs(reordered$A - oracle[perm, , drop = FALSE])), 2e-11)
})

test_that("mixed class panels retain complete probability support and row shape", {
  set.seed(81)
  n <- 36L
  x <- data.frame(x = seq(-1, 1, length.out = n),
    f = factor(rep(letters[1:3], length.out = n)))
  y <- factor(rep(c("b", "a", "c", "a"), length.out = n))
  b <- npcdensbw(xdat = x, ydat = y, bws = c(.15, .5, .2),
    bandwidth.compute = FALSE)
  fit <- npconmode(b, txdat = x, tydat = y, probabilities = TRUE, se = TRUE)
  seed <- .Random.seed
  z <- plot(fit, errors = "asymptotic", band = "all", output = "data",
    view = "fixed", neval = 5, level = "b", gaussian.draws = 9999L)
  expect_identical(.Random.seed, seed)
  expect_length(z, 2L)
  expect_equal(nrow(z[[1]]), 5L)
  expect_equal(nrow(z[[2]]), 3L)
  expect_true(all(vapply(z, function(v) all(is.finite(v$simultaneous.lower)), logical(1))))
  marginal <- plot(fit, errors = "asymptotic", band = "pointwise",
    output = "data", view = "fixed", neval = 5, level = "b")
  for (j in seq_along(z)) expect_equal(z[[j]]$pointwise.lower, marginal[[j]]$lower)
})

test_that("copula ratio factor retains numerator and marginal covariance", {
  set.seed(11)
  d <- data.frame(x = rnorm(45), y = rnorm(45))
  b <- npudensbw(dat = d, bws = c(.5, .6), bandwidth.compute = FALSE)
  obj <- npcopula(b, data = d, u = matrix(rep(c(.3, .5, .7), 2), ncol = 2),
    n.quasi.inv = 30, se = TRUE)
  ex <- .npcopula_eval_xgrid(obj)
  K1 <- outer(ex$x, d$x, function(a, b) dnorm((a - b) / .5) / .5)
  K2 <- outer(ex$y, d$y, function(a, b) dnorm((a - b) / .6) / .6)
  m1 <- rowMeans(K1); m2 <- rowMeans(K2); N <- rowMeans(K1 * K2)
  phi <- (K1 * K2 - N) / (m1 * m2) -
    N / (m1 * m2) * ((K1 - m1) / m1 + (K2 - m2) / m2)
  A <- (phi - rowMeans(phi)) / nrow(d)
  actual <- .np_plot_copula_influence(obj, d, ex, se(obj))
  expect_equal(tcrossprod(actual$A), tcrossprod(A), tolerance = 2e-11)
  expect_equal(sqrt(rowSums(A^2)), se(obj), tolerance = 2e-11)
})
