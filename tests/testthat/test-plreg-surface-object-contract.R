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

test_that("native PL plot reconstruction preserves non-syntactic names", {
  set.seed(892)
  x <- data.frame(value = rnorm(63)); z <- data.frame(value = runif(63))
  names(x) <- "linear x"; names(z) <- "smooth z"
  y <- 2*x[[1]] + cos(z[[1]]) + rnorm(63, sd = .2)
  b <- npplregbw(xdat = x, zdat = z, ydat = y,
                 bws = matrix(.6, 2, 1), bandwidth.compute = FALSE)
  fit <- npplreg(b, se = TRUE)
  for (obj in list(b, fit)) for (surface in c(FALSE, TRUE)) {
    implicit <- plot(obj, plot.behavior = "data", perspective = surface,
                     neval = 5L, errors = "none")
    explicit <- plot(b, xdat = x, ydat = y, zdat = z,
                     plot.behavior = "data", perspective = surface,
                     neval = 5L, errors = "none")
    for (j in seq_along(implicit)) {
      one <- implicit[[j]]; ref <- explicit[[j]]
      expect_identical(names(one$evalx), names(x))
      expect_identical(names(one$evalz), names(z))
      expect_identical(names(coef(one)), names(x))
      expect_equal(unname(coef(one)["linear x"]), unname(coef(fit)), tolerance = 0)
      expect_equal(one$evalx, ref$evalx, tolerance = 0)
      expect_equal(one$evalz, ref$evalz, tolerance = 0)
      expect_equal(fitted(one), fitted(ref), tolerance = 0)
      expect_output(summary(one), "Partially Linear Model")
    }
  }
  names(x) <- "replacement x"; names(z) <- "replacement z"
  replacement <- plot(b, xdat = x, ydat = y, zdat = z,
    plot.behavior = "data", perspective = TRUE, neval = 5L, errors = "none")$r1
  expect_identical(names(replacement$evalx), names(x))
  expect_identical(names(replacement$evalz), names(z))
})
