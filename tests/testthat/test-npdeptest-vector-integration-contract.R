test_that("npdeptest vector integration preserves a skewed dependence oracle", {
  skip_if_not(isTRUE(getOption("npRmpi.pool.active", FALSE)))
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)

  set.seed(84005)
  x <- rexp(20) - 1
  y <- 0.6 * x + 0.4 * rt(20, df = 5)
  expected.bootstrap <- c(
    0.00853060225465798, 0.00397987443375825, 0.00762873928909601,
    0.0121667958169233, 0.00607623631257945, 0.0108004358878692,
    0.00628193347048892, 0.00306063732243844, 0.0093193651817204
  )

  out <- npdeptest(
    data.x = x,
    data.y = y,
    method = "integration",
    bootstrap = TRUE,
    boot.num = 9,
    random.seed = 184005
  )

  expect_lte(abs(out$Srho - 0.0577024275025102), 1e-9)
  expect_lte(max(abs(out$Srho.bootstrap.vec - expected.bootstrap)), 1e-9)
  expect_identical(out$P, 0)
  expect_equal(out$bw.data.x, 0.95821067417595, tolerance = 1e-14)
  expect_equal(out$bw.data.y, 0.476219867158522, tolerance = 1e-14)
  expect_equal(
    out$bw.joint,
    c(0.856724238900334, 0.44585216535079),
    tolerance = 1e-14
  )
})
