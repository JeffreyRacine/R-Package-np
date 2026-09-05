test_that("npsdeptest vector integration preserves a skewed-series oracle", {
  skip_if_not(isTRUE(getOption("npRmpi.pool.active", FALSE)))
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)

  # The oracle also fixes the independent C bandwidth-search seed.
  npseed(42)
  set.seed(84004)
  x <- rexp(20) - 1 + 0.25 * rt(20, df = 5)
  expected.bootstrap <- c(
    0.00765792704386238, 0.00153415516298009, 0.00643904665444113,
    0.00742057756031198, 0.00336669996687301, 0.00284925725849309,
    0.00634086244940879, 0.00169753885446266, 0.00198952590559565
  )

  out <- npsdeptest(
    data = x, lag.num = 1, method = "integration", bootstrap = TRUE,
    B = 9, random.seed = 184004
  )

  expect_lte(abs(out$Srho - 0.00563360665739774), 1e-10)
  expect_lte(max(abs(as.numeric(out$Srho.bootstrap.mat) -
                     expected.bootstrap)), 1e-10)
  expect_identical(dim(out$Srho.bootstrap.mat), c(9L, 1L))
  expect_identical(out$P, 4 / 9)
  expect_identical(out$P.cumulant, 4 / 9)
  # The statistic and bootstrap matrix above are the numerical oracle.
  # Bandwidth endpoints use a tight optimizer-scale guide because harmless
  # floating-point reordering can move the optimizer by a few ulps.
  expect_equal(out$bw.y, 0.802635304793141, tolerance = 5e-11)
  expect_equal(out$bw.y.lag, 0.809899270084268, tolerance = 5e-11)
})
