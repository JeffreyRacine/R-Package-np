# Scalar conditional CVLS expectations use independent dbeta/pbeta weights
# on explicitly deleted samples, including self-dominant boundary rows.
test_that("scalar beta CDF rows retain donors when self weight dominates", {
  old <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old), add = TRUE)
  set.seed(4)
  x <- c(.3, seq(.03, .95, length.out = 19L))
  y <- pmin(pmax(.3 + .4*x + rnorm(20, sd = .1), .02), .98)
  X <- data.frame(x = x); Y <- data.frame(y = y)
  grid <- c(.2, .4, .6, .8)
  for (hx in c(.03, .02, .01)) {
    b <- npcdistbw(xdat = X, ydat = Y, bws = c(.1, hx),
      bandwidth.compute = FALSE, regtype = "lc",
      cxkertype = "beta", cykertype = "beta",
      cxkerbound = "fixed", cxkerlb = 0, cxkerub = 1,
      cykerbound = "fixed", cykerlb = 0, cykerub = 1)
    reference <- mean(vapply(seq_along(x), function(i) {
      logw <- dbeta(x[-i], 1+x[i]/hx^2, 1+(1-x[i])/hx^2, log = TRUE)
      w <- exp(logw - max(logw)); w <- w/sum(w)
      fy <- vapply(grid, function(g)
        sum(w * pbeta(g, 1+y[-i]/.1^2, 1+(1-y[-i])/.1^2)), 0)
      mean(((y[i] <= grid) - fy)^2)
    }, 0))
    observed <- npRmpi:::.npcdistbw_eval_only(
      X, Y, bws = b, gydat = data.frame(y = grid),
      invalid.penalty = "dbmax")$objective
    expect_equal(observed, reference, tolerance = 1e-10)
  }
  hx <- .02; hy <- .1
  b <- npcdensbw(xdat = X, ydat = Y, bws = c(hy, hx),
    bandwidth.compute = FALSE, regtype = "lc", bwmethod = "cv.ls",
    cxkertype = "beta", cykertype = "beta",
    cxkerbound = "fixed", cxkerlb = 0, cxkerub = 1,
    cykerbound = "fixed", cykerlb = 0, cykerub = 1,
    cvls.quadrature.grid = "uniform", cvls.quadrature.points = c(31L, 11L))
  nodes <- seq(0, 1, length.out = 31L)
  quadrature <- rep(1/30, 31L); quadrature[c(1L, 31L)] <- 1/60
  # Beta PDF is query-centred; the CDF above is observation-centred.
  reference <- mean(vapply(seq_along(x), function(i) {
    logw <- dbeta(x[-i], 1+x[i]/hx^2, 1+(1-x[i])/hx^2, log = TRUE)
    w <- exp(logw - max(logw)); w <- w/sum(w)
    fy <- vapply(c(y[i], nodes), function(g)
      sum(w * dbeta(y[-i], 1+g/hy^2, 1+(1-g)/hy^2)), 0)
    2*fy[1L] - sum(quadrature * fy[-1L]^2)
  }, 0))
  observed <- npRmpi:::.npcdensbw_eval_only(
    X, Y, bws = b, invalid.penalty = "dbmax")$objective
  expect_equal(observed, reference, tolerance = 1e-10)

})

# CF-232: evaluate the existing LP solve on omitted weights. The tolerances
# below are the maintainer-approved residuals for these ill-conditioned rows;
# they do not change the production solve guard or other tests' tolerances.
test_that("beta local-linear CVLS uses directly deleted influence rows", {
  old <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old), add = TRUE)
  for (type in c("fixed", "generalized_nn")) {
    if (type == "fixed") {
      set.seed(4)
      x <- c(.3, seq(.03, .95, length.out = 19L))
      y <- pmin(pmax(.3 + .4*x + rnorm(20, sd = .1), .02), .98)
      bw <- c(.1, .02)
      grid <- c(.2, .4, .6, .8)
      tolerance <- 1e-7
    } else {
      set.seed(3)
      x <- rbeta(80, .5, 2)
      y <- pmin(pmax(.3 + .4*sqrt(x) + rnorm(80, sd = .1), .01), .99)
      bw <- c(10, 8)
      grid <- as.numeric(quantile(y, c(.1, .3, .5, .7, .9)))
      tolerance <- 1e-8
    }
    X <- data.frame(x = x); Y <- data.frame(y = y)
    G <- data.frame(y = grid)
    b <- npcdistbw(xdat = X, ydat = Y, bws = bw, bwtype = type,
      bandwidth.compute = FALSE, regtype = "ll",
      cxkertype = "beta", cykertype = "beta",
      cxkerbound = "fixed", cxkerlb = 0, cxkerub = 1,
      cykerbound = "fixed", cykerlb = 0, cykerub = 1)
    reference <- mean(vapply(seq_along(x), function(i) {
      pred <- fitted(npcdist(bws = b, txdat = X[-i, , drop = FALSE],
        tydat = Y[-i, , drop = FALSE],
        exdat = data.frame(x = rep(x[i], length(grid))), eydat = G))
      mean(((y[i] <= grid) - pred)^2)
    }, 0))
    observed <- npRmpi:::.npcdistbw_eval_only(X, Y, bws = b, gydat = G,
      invalid.penalty = "dbmax")$objective
    expect_true(is.finite(observed))
    expect_lt(abs(observed - reference), tolerance)
  }
})
