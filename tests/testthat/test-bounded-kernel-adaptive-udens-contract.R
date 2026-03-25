library(np)

test_that("bounded adaptive_nn is available for unconditional density", {
  set.seed(20260325)
  x <- runif(36)
  xy <- data.frame(x = x)

  bw.fixed <- npudensbw(
    dat = xy,
    bwmethod = "cv.ml",
    bwtype = "fixed",
    ckerbound = "range",
    nmulti = 1
  )
  bw.generalized <- npudensbw(
    dat = xy,
    bwmethod = "cv.ml",
    bwtype = "generalized_nn",
    ckerbound = "range",
    nmulti = 1
  )
  bw.adaptive <- npudensbw(
    dat = xy,
    bwmethod = "cv.ml",
    bwtype = "adaptive_nn",
    ckerbound = "range",
    nmulti = 1
  )

  fit.adaptive <- npudens(bws = bw.adaptive, tdat = xy)

  expect_true(all(is.finite(as.numeric(bw.adaptive$bw))))
  expect_true(is.finite(bw.adaptive$fval))
  expect_true(all(is.finite(as.numeric(fit.adaptive$dens))))

  expect_true(all(is.finite(as.numeric(bw.fixed$bw))))
  expect_true(is.finite(bw.fixed$fval))
  expect_true(all(is.finite(as.numeric(bw.generalized$bw))))
  expect_true(is.finite(bw.generalized$fval))

  expect_false(isTRUE(all.equal(as.numeric(bw.fixed$bw), as.numeric(bw.adaptive$bw))))
  expect_false(isTRUE(all.equal(as.numeric(bw.generalized$bw), as.numeric(bw.adaptive$bw))))
})
