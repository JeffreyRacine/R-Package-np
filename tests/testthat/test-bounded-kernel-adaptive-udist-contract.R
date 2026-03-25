test_that("bounded adaptive_nn is available for unconditional distribution", {
  set.seed(20260325)
  n <- 40
  x <- runif(n)
  xdf <- data.frame(x = x)

  bw.fixed <- npudistbw(
    dat = xdf,
    bwmethod = "cv.cdf",
    bwtype = "fixed",
    ckerorder = 2,
    ckerbound = "range",
    nmulti = 1,
    remin = FALSE,
    itmax = 50,
    tol = 0.1
  )

  bw.gnn <- npudistbw(
    dat = xdf,
    bwmethod = "cv.cdf",
    bwtype = "generalized_nn",
    ckerorder = 2,
    ckerbound = "range",
    nmulti = 1,
    remin = FALSE,
    itmax = 50,
    tol = 0.1
  )

  bw.ad <- npudistbw(
    dat = xdf,
    bwmethod = "cv.cdf",
    bwtype = "adaptive_nn",
    ckerorder = 2,
    ckerbound = "range",
    nmulti = 1,
    remin = FALSE,
    itmax = 50,
    tol = 0.1
  )

  fit.ad <- npudist(bws = bw.ad)
  eval.ad <- predict(fit.ad, newdata = data.frame(x = sort(x)))

  expect_true(is.finite(bw.fixed$fval))
  expect_true(is.finite(bw.gnn$fval))
  expect_true(is.finite(bw.ad$fval))
  expect_true(all(is.finite(as.numeric(bw.ad$bw))))
  expect_true(all(is.finite(as.numeric(fitted(fit.ad)))))
  expect_true(all(is.finite(as.numeric(eval.ad))))
  expect_true(all(as.numeric(eval.ad) >= -0.05))
  expect_true(all(as.numeric(eval.ad) <= 1.05))
  expect_true(all(diff(as.numeric(eval.ad)) >= -1e-8))
  expect_false(isTRUE(all.equal(as.numeric(bw.ad$bw), as.numeric(bw.fixed$bw))))
})
