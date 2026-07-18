test_that("npsdeptest vector integration preserves a skewed-series oracle", {
  skip_if_not(isTRUE(getOption("npRmpi.pool.active", FALSE)))
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)

  set.seed(84004)
  x <- rexp(20) - 1 + 0.25 * rt(20, df = 5)
  expected.bootstrap <- c(
    0.00765793447646022, 0.00153415714572756, 0.00643904204910143,
    0.00742057868607986, 0.00336669954750157, 0.00284925889569732,
    0.00634086050185455, 0.00169752807402632, 0.00198953097270689
  )

  out <- npsdeptest(
    data = x, lag.num = 1, method = "integration", bootstrap = TRUE,
    boot.num = 9, random.seed = 184004
  )

  expect_lte(abs(out$Srho - 0.00563360528836586), 1e-9)
  expect_lte(max(abs(as.numeric(out$Srho.bootstrap.mat) -
                     expected.bootstrap)), 1e-9)
  expect_identical(dim(out$Srho.bootstrap.mat), c(9L, 1L))
  expect_identical(out$P, 4 / 9)
  expect_identical(out$P.cumulant, 4 / 9)
  expect_equal(out$bw.y, 0.802635304793141, tolerance = 1e-14)
  expect_equal(out$bw.y.lag, 0.809899270084268, tolerance = 1e-14)
})
