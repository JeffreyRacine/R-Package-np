test_that("bounded adaptive_nn is available for conditional density", {
  set.seed(20260325)
  n <- 36
  x <- runif(n)
  y <- runif(n)
  xdf <- data.frame(x = x)
  ydf <- data.frame(y = y)

  bw.fixed <- npcdensbw(
    xdat = xdf,
    ydat = ydf,
    bwmethod = "cv.ml",
    bwtype = "fixed",
    cxkerorder = 2,
    cykerorder = 2,
    cxkerbound = "range",
    cykerbound = "range",
    nmulti = 1,
    powell.remin = FALSE,
    itmax = 50,
    tol = 0.1
  )

  bw.gnn <- npcdensbw(
    xdat = xdf,
    ydat = ydf,
    bwmethod = "cv.ml",
    bwtype = "generalized_nn",
    cxkerorder = 2,
    cykerorder = 2,
    cxkerbound = "range",
    cykerbound = "range",
    nmulti = 1,
    powell.remin = FALSE,
    itmax = 50,
    tol = 0.1
  )

  bw.ad <- npcdensbw(
    xdat = xdf,
    ydat = ydf,
    bwmethod = "cv.ml",
    bwtype = "adaptive_nn",
    cxkerorder = 2,
    cykerorder = 2,
    cxkerbound = "range",
    cykerbound = "range",
    nmulti = 1,
    powell.remin = FALSE,
    itmax = 50,
    tol = 0.1
  )

  fit.ad <- npcdens(bws = bw.ad)

  expect_true(is.finite(bw.fixed$fval))
  expect_true(is.finite(bw.gnn$fval))
  expect_true(is.finite(bw.ad$fval))
  expect_true(all(is.finite(as.numeric(bw.ad$xbw))))
  expect_true(all(is.finite(as.numeric(bw.ad$ybw))))
  expect_true(all(is.finite(as.numeric(fitted(fit.ad)))))
  expect_true(all(as.numeric(fitted(fit.ad)) >= 0))
  expect_false(isTRUE(all.equal(as.numeric(bw.ad$xbw), as.numeric(bw.fixed$xbw))))
})
