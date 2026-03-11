test_that("plot contract: bootstrap args require explicit bootstrap mode across engine families (npRmpi)", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  old.auto <- getOption("npRmpi.autodispatch", FALSE)
  on.exit(options(npRmpi.autodispatch = old.auto), add = TRUE)
  on.exit(close_mpi_slaves(), add = TRUE)
  options(npRmpi.autodispatch = TRUE)

  msg <- "plot.errors.method must be set to 'bootstrap' when bootstrap error arguments are supplied"

  set.seed(812)
  n <- 50
  x <- runif(n)
  x2 <- runif(n)
  z <- runif(n)
  y <- sin(2 * pi * x) + 0.3 * z + rnorm(n, sd = 0.08)

  xdat <- data.frame(x = x)
  ydat <- data.frame(y = y)
  idxdat <- data.frame(y = y, x = x, x2 = x2)

  rbw <- npregbw(xdat = xdat, ydat = y, bws = 0.25, bandwidth.compute = FALSE)
  ubw <- npudensbw(dat = xdat, bws = 0.25, bandwidth.compute = FALSE)
  dbw <- npudistbw(dat = xdat, bws = 0.25, bandwidth.compute = FALSE)
  cbw <- npcdensbw(xdat = xdat, ydat = ydat, bws = c(0.30, 0.30), bandwidth.compute = FALSE)
  cdbw <- npcdistbw(xdat = xdat, ydat = ydat, bws = c(0.30, 0.30), bandwidth.compute = FALSE)
  pbw <- npplregbw(
    xdat = data.frame(x = z),
    zdat = data.frame(z = x),
    ydat = y,
    bws = matrix(c(0.30, 0.30), nrow = 2),
    bandwidth.compute = FALSE
  )
  sbw <- npindexbw(
    y ~ x + x2,
    data = idxdat,
    bws = c(1, 0.30, 0.30),
    bandwidth.compute = FALSE,
    method = "ichimura"
  )
  scbw <- npscoefbw(
    xdat = data.frame(x = x),
    zdat = data.frame(z = z),
    ydat = y,
    bws = 0.20,
    bandwidth.compute = FALSE
  )

  expect_error(
    suppressWarnings(plot(rbw, plot.behavior = "data", perspective = FALSE, plot.errors.boot.num = 5)),
    msg,
    fixed = TRUE
  )
  expect_error(
    suppressWarnings(plot(ubw, plot.behavior = "data", perspective = FALSE, plot.errors.boot.num = 5)),
    msg,
    fixed = TRUE
  )
  expect_error(
    suppressWarnings(plot(dbw, plot.behavior = "data", perspective = FALSE, plot.errors.boot.num = 5)),
    msg,
    fixed = TRUE
  )
  expect_error(
    suppressWarnings(plot(cbw, plot.behavior = "data", perspective = FALSE, plot.errors.boot.num = 5)),
    msg,
    fixed = TRUE
  )
  expect_error(
    suppressWarnings(plot(cdbw, plot.behavior = "data", perspective = FALSE, plot.errors.boot.num = 5)),
    msg,
    fixed = TRUE
  )
  expect_error(
    suppressWarnings(
      plot(
        pbw,
        xdat = data.frame(x = z),
        ydat = y,
        zdat = data.frame(z = x),
        plot.behavior = "data",
        perspective = FALSE,
        plot.errors.boot.num = 5
      )
    ),
    msg,
    fixed = TRUE
  )
  expect_error(
    suppressWarnings(
      plot(
        sbw,
        xdat = data.frame(x = x, x2 = x2),
        ydat = y,
        plot.behavior = "data",
        perspective = FALSE,
        plot.errors.boot.num = 5
      )
    ),
    msg,
    fixed = TRUE
  )
  expect_error(
    suppressWarnings(
      plot(
        scbw,
        xdat = data.frame(x = x),
        ydat = y,
        zdat = data.frame(z = z),
        plot.behavior = "data",
        perspective = FALSE,
        plot.errors.boot.num = 5
      )
    ),
    msg,
    fixed = TRUE
  )
})
