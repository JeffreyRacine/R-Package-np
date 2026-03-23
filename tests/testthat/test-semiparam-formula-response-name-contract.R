test_that("semiparametric formula fits preserve response-name metadata for plotting", {
  set.seed(123)
  n <- 50L

  dat_pl <- data.frame(x = rnorm(n), z = sort(runif(n)))
  dat_pl$y <- 1 + dat_pl$x + sin(2 * pi * dat_pl$z) + rnorm(n, sd = 0.08)
  fit_pl <- np::npplreg(
    y ~ x | z,
    data = dat_pl,
    nmulti = 1,
    regtype = "lp",
    degree = 1,
    bwtype = "adaptive_nn"
  )
  expect_identical(fit_pl$bws$ynames, "y")

  dat_sc <- data.frame(x = runif(n), z = sort(runif(n)))
  dat_sc$y <- (1 + dat_sc$z^2) * dat_sc$x + rnorm(n, sd = 0.08)
  fit_sc <- np::npscoef(
    y ~ x | z,
    data = dat_sc,
    nmulti = 1,
    regtype = "lp",
    degree = 1,
    bwtype = "adaptive_nn",
    errors = FALSE,
    betas = FALSE
  )
  expect_identical(fit_sc$bws$ynames, "y")

  dat_si <- data.frame(x1 = runif(n, -1, 1), x2 = runif(n, -1, 1))
  idx <- dat_si$x1 + 0.5 * dat_si$x2
  dat_si$y <- sin(idx) + 0.25 * idx^2 + rnorm(n, sd = 0.05)
  fit_si <- np::npindex(
    y ~ x1 + x2,
    data = dat_si,
    method = "ichimura",
    nmulti = 1,
    regtype = "lp",
    degree = 1,
    bwtype = "adaptive_nn"
  )
  expect_identical(fit_si$bws$ynames, "y")
})

test_that("semiparametric formula-fit plots do not regress on response labels", {
  set.seed(456)
  n <- 40L

  dat_pl <- data.frame(x = rnorm(n), z = sort(runif(n)))
  dat_pl$y <- 1 + dat_pl$x + sin(2 * pi * dat_pl$z) + rnorm(n, sd = 0.08)
  fit_pl <- np::npplreg(
    y ~ x | z,
    data = dat_pl,
    nmulti = 1,
    regtype = "lp",
    degree = 1,
    bwtype = "adaptive_nn"
  )

  dat_sc <- data.frame(x = runif(n), z = sort(runif(n)))
  dat_sc$y <- (1 + dat_sc$z^2) * dat_sc$x + rnorm(n, sd = 0.08)
  fit_sc <- np::npscoef(
    y ~ x | z,
    data = dat_sc,
    nmulti = 1,
    regtype = "lp",
    degree = 1,
    bwtype = "adaptive_nn",
    errors = FALSE,
    betas = FALSE
  )

  dat_si <- data.frame(x1 = runif(n, -1, 1), x2 = runif(n, -1, 1))
  idx <- dat_si$x1 + 0.5 * dat_si$x2
  dat_si$y <- sin(idx) + 0.25 * idx^2 + rnorm(n, sd = 0.05)
  fit_si <- np::npindex(
    y ~ x1 + x2,
    data = dat_si,
    method = "ichimura",
    nmulti = 1,
    regtype = "lp",
    degree = 1,
    bwtype = "adaptive_nn"
  )

  path <- tempfile(fileext = ".png")
  grDevices::png(path, width = 800, height = 600)
  on.exit({
    grDevices::dev.off()
    unlink(path)
  }, add = TRUE)

  expect_no_error(plot(fit_pl))
  expect_no_error(plot(fit_sc))
  expect_no_error(plot(fit_si))
})
