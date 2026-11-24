test_that("pooled quadratic Gaussian slope hats retain the fitted ridge", {
  old <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old), add = TRUE)
  # At radius one, only x=0 and x=1 have nonzero Gaussian weights.
  # Their raw linear and quadratic columns coincide exactly.
  x <- data.frame(x = c(0, 0, 1, 60:84),
                  z = factor(rep(c("a", "b"), 14L)),
                  w = factor(rep(c("c", "c", "d", "d"), 7L)))
  y <- sin(x$x)
  ev <- x[1L, , drop = FALSE]
  for (k in c(3L, 21L)) {
    b <- npregbw(xdat = x, ydat = y, bws = c(k, .5, .5),
      bwtype = "generalized_nn", ckertype = "gaussian", ckerbound = "none",
      ukertype = "aitchisonaitken", regtype = "lp", basis = "glp", degree = 2L,
      bernstein.basis = FALSE, bandwidth.compute = FALSE)
    H <- npreghat(b, txdat = x, exdat = ev, s = 1L)
    ridge <- attr(H, "ridge.used", exact = TRUE)
    expect_identical(ridge > 0, k == 3L)
    fit <- npreg(b, txdat = x, tydat = y, exdat = ev, gradients = TRUE, se = FALSE)
    expect_equal(as.numeric(H %*% y), unname(gradients(fit)[1L, 1L]), tolerance = 1e-8)
    if (k == 3L) {
      a <- npRmpi:::.np_bias_regression_prepare(b, x, y, gradient.order = 1L)
      expect_error(npRmpi:::.np_bias_regression_rows(a, ev, data.frame(row = 1L, target = "d1")),
        "pooled mixed-data correction requires unregularized main rows", fixed = TRUE)
    }
  }
})
