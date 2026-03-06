library(np)

test_that("npscoefbw fast counter fires only on the z-irrelevant sentinel", {
  old_opts <- options(
    np.messages = FALSE,
    np.largeh.rel.tol = 0.05,
    np.disc.upper.rel.tol = 0.05
  )
  on.exit(options(old_opts), add = TRUE)

  set.seed(77)
  n <- 200L
  x <- runif(n)
  z <- runif(n)

  dat_lc <- data.frame(
    y = rnorm(n, sd = 0.5 * sd(x)),
    x = x,
    z = z
  )
  dat_ll <- data.frame(
    y = x + rnorm(n, sd = 0.5 * sd(x)),
    x = x,
    z = z
  )

  bw_lc <- npscoefbw(y ~ x | z, data = dat_lc, regtype = "lc", bwmethod = "cv.ls", nmulti = 1)
  bw_ll <- npscoefbw(y ~ x | z, data = dat_ll, regtype = "ll", bwmethod = "cv.ls", nmulti = 1)
  bw_lp <- npscoefbw(
    y ~ x | z,
    data = dat_ll,
    regtype = "lp",
    degree = 1,
    bwmethod = "cv.ls",
    nmulti = 1
  )

  expect_true(is.finite(as.numeric(bw_lc$num.feval.fast)))
  expect_true(is.finite(as.numeric(bw_ll$num.feval.fast)))
  expect_true(is.finite(as.numeric(bw_lp$num.feval.fast)))

  expect_gt(as.integer(bw_lc$num.feval.fast), 0L)
  expect_gt(as.integer(bw_ll$num.feval.fast), 0L)
  expect_gt(as.integer(bw_lp$num.feval.fast), 0L)

  expect_equal(as.numeric(bw_ll$fval), as.numeric(bw_lp$fval), tolerance = 1e-10)
  expect_equal(as.integer(bw_ll$num.feval.fast), as.integer(bw_lp$num.feval.fast))
})
