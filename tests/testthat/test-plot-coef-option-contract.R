test_that("plot coef option: npscoef supports coef=TRUE in data mode", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")

  old.auto <- getOption("npRmpi.autodispatch", FALSE)
  on.exit(options(npRmpi.autodispatch = old.auto), add = TRUE)
  options(npRmpi.autodispatch = TRUE)

  set.seed(105)
  n <- 60
  x <- runif(n)
  z <- runif(n, -2, 2)
  y <- x * exp(z) * (1 + rnorm(n, sd = 0.15))

  fit <- npscoef(y ~ x | z, regtype = "ll", betas = TRUE)
  out <- suppressWarnings(
    plot(
      fit,
      coef = TRUE,
      coef.index = 1,
      perspective = FALSE,
      neval = 20,
      plot.behavior = "plot-data",
      plot.errors.method = "none"
    )
  )

  expect_type(out, "list")
  expect_true(length(out) > 0)
  expect_true(all(vapply(out, inherits, logical(1), "smoothcoefficient")))
})

test_that("plot coef option: npplreg supports coef=TRUE plot-data payload", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")

  old.auto <- getOption("npRmpi.autodispatch", FALSE)
  on.exit(options(npRmpi.autodispatch = old.auto), add = TRUE)
  options(npRmpi.autodispatch = TRUE)

  set.seed(106)
  n <- 80
  x <- runif(n)
  z <- runif(n, -2, 2)
  y <- 1 + 0.7 * x + sin(z) + rnorm(n, sd = 0.15)
  xdat <- data.frame(x = x)
  zdat <- data.frame(z = z)

  fit <- npplreg(y ~ x | z, regtype = "ll")
  out <- suppressWarnings(
    plot(
      fit,
      xdat = xdat,
      ydat = y,
      zdat = zdat,
      coef = TRUE,
      plot.behavior = "plot-data",
      plot.errors.method = "none"
    )
  )

  expect_type(out, "list")
  expect_true(all(c("coefficients", "coefficient.stderr", "fit") %in% names(out)))
  expect_true(is.numeric(out$coefficients))
  expect_true(length(out$coefficients) >= 1L)
})
