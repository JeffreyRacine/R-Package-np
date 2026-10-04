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
