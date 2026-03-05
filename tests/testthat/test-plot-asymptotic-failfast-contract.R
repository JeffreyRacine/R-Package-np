test_that("plot contract: asymptotic mode fails fast for plbandwidth and sibandwidth (npRmpi)", {
  skip_if_not_installed("np")
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  old.auto <- getOption("npRmpi.autodispatch", FALSE)
  on.exit(options(npRmpi.autodispatch = old.auto), add = TRUE)
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  options(npRmpi.autodispatch = TRUE)

  set.seed(922)
  n <- 65
  x <- runif(n)
  z <- runif(n)
  x2 <- runif(n)
  y <- sin(2 * pi * z) + 0.4 * x + rnorm(n, sd = 0.08)

  pbw <- np::npplregbw(
    xdat = data.frame(x = x),
    zdat = data.frame(z = z),
    ydat = y,
    bws = matrix(c(0.30, 0.30), nrow = 2),
    bandwidth.compute = FALSE
  )
  sbw <- np::npindexbw(
    y ~ x + x2,
    data = data.frame(y = y, x = x, x2 = x2),
    bws = c(1, 0.30, 0.30),
    bandwidth.compute = FALSE,
    method = "ichimura"
  )

  expect_error(
    suppressWarnings(
      plot(
        pbw,
        xdat = data.frame(x = x),
        ydat = y,
        zdat = data.frame(z = z),
        plot.behavior = "data",
        plot.errors.method = "asymptotic",
        perspective = FALSE
      )
    ),
    "asymptotic errors are unsupported for partially linear regression plots",
    fixed = TRUE
  )

  expect_error(
    suppressWarnings(
      plot(
        sbw,
        xdat = data.frame(x = x, x2 = x2),
        ydat = y,
        plot.behavior = "data",
        plot.errors.method = "asymptotic",
        perspective = FALSE
      )
    ),
    "asymptotic errors are unsupported for single-index regression plots",
    fixed = TRUE
  )
})
