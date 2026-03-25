test_that("bounded adaptive_nn is available for conditional distribution", {
  set.seed(20260325)
  n <- 36
  x <- runif(n)
  y <- runif(n)
  xdf <- data.frame(x = x)
  ydf <- data.frame(y = y)

  bw.fixed <- npcdistbw(
    xdat = xdf,
    ydat = ydf,
    bwmethod = "cv.ls",
    bwtype = "fixed",
    cxkerorder = 2,
    cykerorder = 2,
    cxkerbound = "range",
    cykerbound = "range",
    nmulti = 1,
    remin = FALSE,
    itmax = 50,
    tol = 0.1
  )

  bw.gnn <- npcdistbw(
    xdat = xdf,
    ydat = ydf,
    bwmethod = "cv.ls",
    bwtype = "generalized_nn",
    cxkerorder = 2,
    cykerorder = 2,
    cxkerbound = "range",
    cykerbound = "range",
    nmulti = 1,
    remin = FALSE,
    itmax = 50,
    tol = 0.1
  )

  bw.ad <- npcdistbw(
    xdat = xdf,
    ydat = ydf,
    bwmethod = "cv.ls",
    bwtype = "adaptive_nn",
    cxkerorder = 2,
    cykerorder = 2,
    cxkerbound = "range",
    cykerbound = "range",
    nmulti = 1,
    remin = FALSE,
    itmax = 50,
    tol = 0.1
  )

  fit.ad <- npcdist(bws = bw.ad)

  expect_true(is.finite(bw.fixed$fval))
  expect_true(is.finite(bw.gnn$fval))
  expect_true(is.finite(bw.ad$fval))
  expect_true(all(is.finite(as.numeric(bw.ad$xbw))))
  expect_true(all(is.finite(as.numeric(bw.ad$ybw))))
  expect_true(all(is.finite(as.numeric(fitted(fit.ad)))))
  expect_true(all(as.numeric(fitted(fit.ad)) >= -0.05))
  expect_true(all(as.numeric(fitted(fit.ad)) <= 1.10))
  expect_false(isTRUE(all.equal(as.numeric(bw.ad$xbw), as.numeric(bw.fixed$xbw))))
})
