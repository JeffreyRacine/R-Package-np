test_that("plot coef option: npscoef supports coef=TRUE in plot-data mode", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  withr::local_pdf(NULL)

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
      output = "plot-data",
      errors = "none"
    )
  )

  expect_type(out, "list")
  expect_true(length(out) > 0)
  expect_s3_class(out, "np_plot_data")
  expect_true(all(out$targets$quantity == "coefficient"))
  expect_true(all(vapply(out$panels, function(p) all(is.finite(p$values$estimate)), logical(1))))
})

test_that("plot coef option: npplreg supports coef=TRUE plot-data payload", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  withr::local_pdf(NULL)

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
      output = "plot-data",
      errors = "none"
    )
  )

  expect_type(out, "list")
  expect_s3_class(out, "np_plot_data")
  expect_identical(out$targets$quantity, "coefficient")
  expect_true(is.numeric(plot_target(out)$estimate))
  expect_gte(nrow(plot_target(out)), 1L)
  expect_equal(plot_target(out)$estimate, unname(coef(fit)))
})
