test_that("npreg formula fits preserve the response name for plotting metadata", {
  set.seed(42)
  n <- 60L
  dat <- data.frame(x = runif(n, -1, 1))
  dat$y <- dat$x^2 + rnorm(n, sd = 0.25 * stats::sd(dat$x))

  fit <- np::npreg(
    y ~ x,
    data = dat,
    nmulti = 1,
    regtype = "lp",
    degree = 1,
    bwtype = "adaptive_nn"
  )

  expect_identical(fit$bws$ynames, "y")

  path <- tempfile(fileext = ".png")
  grDevices::png(path, width = 800, height = 600)
  on.exit({
    grDevices::dev.off()
    unlink(path)
  }, add = TRUE)

  expect_no_error(
    plot(fit, errors = "asymptotic", view = "fixed")
  )
})
