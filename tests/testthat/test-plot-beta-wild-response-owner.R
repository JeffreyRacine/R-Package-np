test_that("beta wild derivatives retain their per-response endpoint calculation", {
  withr::local_options(np.messages = FALSE)
  withr::local_preserve_seed()
  x <- data.frame(x = c(0, .12, .3, .55, .82, 1))
  y <- rep(7, 6L)
  ex <- data.frame(x = c(0, .4, 1))
  reuse <- getFromNamespace(".np_wild_boot_from_reghat_exact", "np")
  testthat::local_mocked_bindings(
    .np_wild_boot_from_reghat_operator = function(...) stop("beta must retain its response owner"),
    .package = "np")
  for (regtype in c("lc", "lp"))
    for (order in c(2L, 4L, 6L, 8L)) {
      b <- npregbw(xdat = x, ydat = y, bws = .16,
        bandwidth.compute = FALSE, ckertype = "beta", ckerorder = order,
        ckerbound = "fixed", ckerlb = 0, ckerub = 1, regtype = regtype,
        degree = if (regtype == "lp") 0L else NULL)
      expected <- gradients(npreg(b, txdat = x, tydat = y,
        exdat = ex, gradients = TRUE))[, 1L]
      set.seed(83)
      invisible(runif(nrow(x) * 7L))
      expected.seed <- .Random.seed
      set.seed(83)
      actual <- reuse(xdat = x, exdat = ex, bws = b, ydat = y,
        fit.mean.train = y, B = 7L, gradients = TRUE)
      expect_identical(actual$t0, as.vector(expected))
      expect_equal(actual$t, matrix(rep(expected, each = 7L), 7L, nrow(ex)),
                   tolerance = 0)
      expect_identical(.Random.seed, expected.seed)
      if (order == 8L) expect_identical(actual$t0, rep(0, 3L))
    }
})
