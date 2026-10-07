test_that("partially linear surface payloads retain named dimensions", {
  set.seed(752)
  n <- 61L
  d <- data.frame(x = rnorm(n), z = runif(n))
  d$y <- 2 * d$x + cos(d$z) + rnorm(n, sd = .2)
  bw <- npplregbw(xdat = d["x"], zdat = d["z"], ydat = d$y,
                  bws = matrix(.55, 2, 1), bandwidth.compute = FALSE)
  fit <- npplreg(bw, se = TRUE)
  for (common in c(FALSE, TRUE)) {
    for (surface in c(FALSE, TRUE)) {
      payload <- plot(fit, plot.behavior = "data", neval = 5L,
                      errors = "none", perspective = surface,
                      common.scale = common)
      for (one in payload) {
        expect_s3_class(one, "plregression")
        expect_identical(one$nobs, if (surface) 25L else 5L)
        expect_identical(one$xndim, 1L)
        expect_identical(one$zndim, 1L)
        expect_identical(one$data.xnames, "x")
        expect_identical(one$data.znames, "z")
        expect_true(all(is.na(unlist(one[c("R2", "MSE", "MAE", "MAPE", "CORR", "SIGN")]))))
        expect_output(print(one), "Partially Linear Model")
        expect_output(summary(one), "Partially Linear Model")
        ref <- npplreg(bw, exdat = one$evalx, ezdat = one$evalz)
        expect_equal(fitted(one), fitted(ref), tolerance = 0)
      }
    }
  }
})
