make_semiparam_formula_fits <- function(seed = 123, n = 24L) {
  set.seed(seed)

  dat_pl <- data.frame(x = rnorm(n), z = sort(runif(n)))
  dat_pl$y <- 1 + dat_pl$x + sin(2 * pi * dat_pl$z) + rnorm(n, sd = 0.08)
  bw_pl <- np::npplregbw(
    y ~ x | z,
    data = dat_pl,
    bws = matrix(c(0.7, 0.5), nrow = 2),
    bandwidth.compute = FALSE,
    regtype = "lp",
    degree = 1,
    bwtype = "fixed"
  )
  fit_pl <- np::npplreg(bws = bw_pl)

  dat_sc <- data.frame(x = runif(n), z = sort(runif(n)))
  dat_sc$y <- (1 + dat_sc$z^2) * dat_sc$x + rnorm(n, sd = 0.08)
  bw_sc <- np::npscoefbw(
    y ~ x | z,
    data = dat_sc,
    bws = 0.5,
    bandwidth.compute = FALSE,
    regtype = "lp",
    degree = 1,
    bwtype = "fixed"
  )
  fit_sc <- np::npscoef(bws = bw_sc, errors = FALSE, betas = FALSE)

  dat_si <- data.frame(x1 = runif(n, -1, 1), x2 = runif(n, -1, 1))
  idx <- dat_si$x1 + 0.5 * dat_si$x2
  dat_si$y <- sin(idx) + 0.25 * idx^2 + rnorm(n, sd = 0.05)
  bw_si <- np::npindexbw(
    y ~ x1 + x2,
    data = dat_si,
    bws = c(1, 0.7, 0.7),
    bandwidth.compute = FALSE,
    method = "ichimura",
    regtype = "lp",
    degree = 1,
    bwtype = "fixed"
  )
  fit_si <- np::npindex(bws = bw_si)

  list(pl = fit_pl, sc = fit_sc, si = fit_si)
}

test_that("semiparametric formula fits preserve response-name metadata for plotting", {
  fits <- make_semiparam_formula_fits(seed = 123)

  expect_identical(fits$pl$bws$ynames, "y")
  expect_identical(fits$sc$bws$ynames, "y")
  expect_identical(fits$si$bws$ynames, "y")
})

test_that("semiparametric formula-fit plots do not regress on response labels", {
  fits <- make_semiparam_formula_fits(seed = 456)
  out.pl <- expect_no_error(plot(fits$pl, plot.behavior = "data", perspective = FALSE))
  out.sc <- expect_no_error(plot(fits$sc, plot.behavior = "data", perspective = FALSE))
  out.si <- expect_no_error(plot(fits$si, plot.behavior = "data", perspective = FALSE))

  expect_type(out.pl, "list")
  expect_type(out.sc, "list")
  expect_type(out.si, "list")
  expect_true(length(out.pl) > 0L)
  expect_true(length(out.sc) > 0L)
  expect_true(length(out.si) > 0L)
})
