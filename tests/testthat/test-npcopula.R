test_that("npcopula basic functionality works", {
  # skip_on_cran()
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")

  data("faithful")
  # npcopula can take a dbandwidth (from npudistbw) for copula 
  # or bandwidth (from npudensbw) for copula density
  
  bw <- npudistbw(dat=faithful, bws=c(0.5, 5), bandwidth.compute=FALSE)
  
  # Copula CDF
  cop <- npcopula(data=faithful, bws=bw)
  expect_s3_class(cop, "npcopula")
  expect_s3_class(cop, "data.frame")
  expect_true("copula" %in% names(cop))
  expect_identical(attr(cop, "target"), "distribution")
  expect_identical(attr(cop, "evaluation"), "sample")
  expect_true(all(cop$copula >= 0 & cop$copula <= 1))
  
  # Copula density
  bw_dens <- npudensbw(dat=faithful, bws=c(0.5, 5), bandwidth.compute=FALSE)
  cop_dens <- npcopula(data=faithful, bws=bw_dens)
  expect_s3_class(cop_dens, "npcopula")
  expect_s3_class(cop_dens, "data.frame")
  expect_true("copula" %in% names(cop_dens))
  expect_identical(attr(cop_dens, "target"), "density")
  expect_identical(attr(cop_dens, "evaluation"), "sample")
  expect_true(all(cop_dens$copula >= 0))
})

test_that("npcopula u-grid path completes with active slaves", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")

  set.seed(1)
  n <- 120L
  mydat <- data.frame(
    x = rnorm(n),
    y = 0.8 * rnorm(n) + rnorm(n, sd = 0.2)
  )
  u <- data.frame(x = c(0.25, 0.5, 0.75), y = c(0.25, 0.5, 0.75))

  bw <- npudistbw(dat = mydat, bws = c(0.8, 0.8), bandwidth.compute = FALSE)
  cop <- npcopula(bws = bw, data = mydat, u = u, n.quasi.inv = 60)

  expect_s3_class(cop, "npcopula")
  expect_s3_class(cop, "data.frame")
  expect_true(all(c("copula", "u1", "u2", "x", "y") %in% names(cop)))
  expect_identical(attr(cop, "target"), "distribution")
  expect_identical(attr(cop, "evaluation"), "grid")
  expect_equal(attr(cop, "grid.dim"), c(3L, 3L))
  expect_equal(nrow(cop), 9L)
  expect_true(all(is.finite(cop$copula)))
})

test_that("npcopula formula route builds a plot-ready grid by default", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")

  set.seed(42)
  dat <- data.frame(x = rnorm(50), y = 0.5 * rnorm(50) + rnorm(50, sd = 0.1))

  fit <- npcopula(~ x + y, data = dat, neval = 4, nmulti = 1)

  expect_s3_class(fit, "npcopula")
  expect_s3_class(fit, "data.frame")
  expect_equal(nrow(fit), 16L)
  expect_identical(attr(fit, "target"), "distribution")
  expect_identical(attr(fit, "evaluation"), "grid")
  expect_true(isTRUE(attr(fit, "u.auto")))
  expect_equal(attr(fit, "grid.dim"), c(4L, 4L))
  expect_true(inherits(attr(fit, "bws"), "dbandwidth"))
})

test_that("npcopula formula route can request sample evaluation and density target", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")

  set.seed(43)
  dat <- data.frame(x = rnorm(40), y = 0.5 * rnorm(40) + rnorm(40, sd = 0.1))

  sample.fit <- npcopula(~ x + y, data = dat, evaluation = "sample", nmulti = 1)
  expect_equal(nrow(sample.fit), nrow(dat))
  expect_identical(attr(sample.fit, "evaluation"), "sample")

  dens.fit <- npcopula(~ x + y, data = dat, target = "density", neval = 3, nmulti = 1)
  expect_equal(nrow(dens.fit), 9L)
  expect_identical(attr(dens.fit, "target"), "density")
  expect_true(inherits(attr(dens.fit, "bws"), "bandwidth"))
})

test_that("npcopula formula route rejects unsafe automatic high-dimensional grids", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")

  dat <- data.frame(x = rnorm(20), y = rnorm(20), z = rnorm(20))

  expect_error(
    npcopula(~ x + y + z, data = dat, nmulti = 1),
    "automatic copula probability grids"
  )
})

test_that("npcopula target conflicts with explicit bandwidth object fail clearly", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")

  data("faithful")
  bw <- npudistbw(dat=faithful, bws=c(0.5, 5), bandwidth.compute=FALSE)
  expect_error(
    npcopula(data = faithful, bws = bw, target = "density"),
    "target.*conflicts"
  )
})

test_that("npcopula fitted and basic plot methods work", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")

  data("faithful")
  bw <- npudistbw(dat=faithful, bws=c(0.5, 5), bandwidth.compute=FALSE)
  u <- data.frame(eruptions = seq(0, 1, length.out = 4),
                  waiting = seq(0, 1, length.out = 4))
  fit <- npcopula(data = faithful, bws = bw, u = u, n.quasi.inv = 40)

  expect_equal(fitted(fit), fit$copula)
  expect_silent(plot(fit))
  expect_silent(plot(fit, col = NULL))
  expect_silent(plot(fit, col = "lightblue", border = "grey40",
                     theta = 20, phi = 25, shade = 0.3))
  expect_silent(plot(fit, view = "image", perspective = FALSE))
  expect_silent(plot(fit, view = "contour", perspective = FALSE))
  expect_true("copula" %in% .np_progress_single_line_surfaces())
  expect_error(plot(fit, view = "image", renderer = "rgl"),
               "view='surface'")

  one.point <- npcopula(data = faithful, bws = bw, u = c(0.5, 0.5), n.quasi.inv = 40)
  expect_equal(nrow(one.point), 1L)
  expect_equal(attr(one.point, "grid.dim"), c(1L, 1L))
})

test_that("npcopula plot intervals work for grid surfaces", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")

  set.seed(44)
  dat <- data.frame(x = rnorm(45), y = rnorm(45))
  fit <- npcopula(~ x + y, data = dat, neval = 3, nmulti = 1)

  asym <- plot(fit, output = "data", errors = "asymptotic", band = "all")
  expect_true(all(c("lower", "upper", "pointwise.lower",
                    "bonferroni.upper") %in% names(asym)))
  expect_equal(nrow(asym), nrow(fit))

  boot <- plot(fit, output = "data", errors = "bootstrap",
               bootstrap = "inid", B = 3, band = "pmzsd")
  expect_true(all(c("center", "lower", "upper") %in% names(boot)))
  expect_equal(nrow(boot), nrow(fit))

  dens.bw <- npudensbw(dat = dat, bws = c(0.7, 0.7),
                       bandwidth.compute = FALSE)
  dens.fit <- npcopula(
    bws = dens.bw,
    data = dat,
    u = data.frame(x = seq(0.2, 0.8, length.out = 3),
                   y = seq(0.2, 0.8, length.out = 3)),
    n.quasi.inv = 40
  )
  dens <- plot(dens.fit, output = "data", errors = "bootstrap",
               bootstrap = "inid", B = 3, band = "pmzsd")
  expect_true(all(c("center", "lower", "upper") %in% names(dens)))
  expect_equal(nrow(dens), nrow(dens.fit))

  sample.fit <- npcopula(~ x + y, data = dat, evaluation = "sample", nmulti = 1)
  expect_error(plot(sample.fit, errors = "asymptotic"), "grid evaluation")
})
