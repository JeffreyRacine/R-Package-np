test_that("bounded adaptive_nn is available for unconditional regression", {
  set.seed(20260325)
  n <- 36
  x <- runif(n)
  y <- sin(2 * pi * x) + rnorm(n, sd = 0.15)
  xdf <- data.frame(x = x)

  bw.fixed <- npregbw(
    xdat = xdf,
    ydat = y,
    regtype = "ll",
    bwmethod = "cv.ls",
    bwtype = "fixed",
    ckerorder = 2,
    ckerbound = "range",
    nmulti = 1,
    remin = FALSE,
    itmax = 50,
    tol = 0.1
  )

  bw.gnn <- npregbw(
    xdat = xdf,
    ydat = y,
    regtype = "ll",
    bwmethod = "cv.ls",
    bwtype = "generalized_nn",
    ckerorder = 2,
    ckerbound = "range",
    nmulti = 1,
    remin = FALSE,
    itmax = 50,
    tol = 0.1
  )

  bw.ad <- npregbw(
    xdat = xdf,
    ydat = y,
    regtype = "ll",
    bwmethod = "cv.ls",
    bwtype = "adaptive_nn",
    ckerorder = 2,
    ckerbound = "range",
    nmulti = 1,
    remin = FALSE,
    itmax = 50,
    tol = 0.1
  )

  fit.ad <- npreg(bws = bw.ad, gradients = TRUE)

  expect_true(is.finite(bw.fixed$fval))
  expect_true(is.finite(bw.gnn$fval))
  expect_true(is.finite(bw.ad$fval))
  expect_true(all(is.finite(as.numeric(bw.ad$bw))))
  expect_true(all(is.finite(as.numeric(fitted(fit.ad)))))
  expect_true(all(is.finite(as.numeric(fit.ad$gradients))))
  expect_false(isTRUE(all.equal(as.numeric(bw.ad$bw), as.numeric(bw.fixed$bw))))
  expect_false(isTRUE(all.equal(as.numeric(bw.ad$bw), as.numeric(bw.gnn$bw))))
})
