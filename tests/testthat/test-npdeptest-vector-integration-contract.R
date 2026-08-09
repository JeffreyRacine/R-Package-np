test_that("npdeptest vector integration preserves a skewed dependence oracle", {
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)

  set.seed(84005)
  x <- rexp(20) - 1
  y <- 0.6 * x + 0.4 * rt(20, df = 5)
  expected.bootstrap <- c(
    0.00855095734386373, 0.00399798973921726, 0.00762873596240952,
    0.0121667954163648, 0.00607623553580284, 0.0108003145195545,
    0.00628197180463641, 0.0030689093983868, 0.00931933980398223
  )

  out <- npdeptest(
    data.x = x,
    data.y = y,
    method = "integration",
    bootstrap = TRUE,
    boot.num = 9,
    random.seed = 184005
  )

  expect_lte(abs(out$Srho - 0.0577023468089846), 1e-10)
  expect_lte(max(abs(out$Srho.bootstrap.vec - expected.bootstrap)), 1e-10)
  expect_identical(out$P, 0)
  expect_equal(out$bw.data.x, 0.95821067417595, tolerance = 1e-14)
  expect_equal(out$bw.data.y, 0.476219867158522, tolerance = 1e-14)
  expect_equal(
    out$bw.joint,
    c(0.856724238900334, 0.44585216535079),
    tolerance = 1e-14
  )
})
