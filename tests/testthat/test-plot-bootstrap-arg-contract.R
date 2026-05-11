test_that("plot contract: bootstrap args require explicit bootstrap mode across engine families (npRmpi)", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  old.auto <- getOption("npRmpi.autodispatch", FALSE)
  on.exit(options(npRmpi.autodispatch = old.auto), add = TRUE)
  on.exit(close_mpi_slaves(), add = TRUE)
  options(npRmpi.autodispatch = TRUE)

  msg <- "bootstrap controls require errors = \"bootstrap\""

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
    suppressWarnings(plot(rbw, output = "data", perspective = FALSE, B = 5)),
    msg,
    fixed = TRUE
  )
  expect_error(
    suppressWarnings(plot(ubw, output = "data", perspective = FALSE, B = 5)),
    msg,
    fixed = TRUE
  )
  expect_error(
    suppressWarnings(plot(dbw, output = "data", perspective = FALSE, B = 5)),
    msg,
    fixed = TRUE
  )
  expect_error(
    suppressWarnings(plot(cbw, output = "data", perspective = FALSE, B = 5)),
    msg,
    fixed = TRUE
  )
  expect_error(
    suppressWarnings(plot(cdbw, output = "data", perspective = FALSE, B = 5)),
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
        output = "data",
        perspective = FALSE,
        B = 5
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
        output = "data",
        perspective = FALSE,
        B = 5
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
        output = "data",
        perspective = FALSE,
        B = 5
      )
    ),
    msg,
    fixed = TRUE
  )
})

test_that("plot contract: package legends are user-controllable on band-all plots (npRmpi)", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  old.auto <- getOption("npRmpi.autodispatch", FALSE)
  on.exit(options(npRmpi.autodispatch = old.auto), add = TRUE)
  on.exit(close_mpi_slaves(), add = TRUE)
  options(npRmpi.autodispatch = TRUE)

  set.seed(20260511)
  n <- 36L
  x <- runif(n)
  y <- sin(2 * pi * x) + rnorm(n, sd = 0.1)
  rbw <- npregbw(
    xdat = data.frame(x = x),
    ydat = y,
    bws = 0.35,
    bandwidth.compute = FALSE
  )

  old.dev <- grDevices::dev.cur()
  grDevices::pdf(file = tempfile(fileext = ".pdf"))
  on.exit({
    grDevices::dev.off()
    if (old.dev > 1L)
      grDevices::dev.set(old.dev)
  }, add = TRUE)

  expect_silent(suppressWarnings(plot(
    rbw, xdat = data.frame(x = x), ydat = y, perspective = FALSE,
    errors = "bootstrap", bootstrap = "wild", B = 9L, band = "all",
    neval = 6L, legend = FALSE
  )))
  expect_silent(suppressWarnings(plot(
    rbw, xdat = data.frame(x = x), ydat = y, perspective = FALSE,
    errors = "bootstrap", bootstrap = "wild", B = 9L, band = "all",
    neval = 6L, legend = "bottomright"
  )))
  expect_silent(suppressWarnings(plot(
    rbw, xdat = data.frame(x = x), ydat = y, perspective = FALSE,
    errors = "bootstrap", bootstrap = "wild", B = 9L, band = "all",
    neval = 6L, legend = list(x = "bottomleft", bty = "o")
  )))
  expect_error(
    suppressWarnings(plot(
      rbw, xdat = data.frame(x = x), ydat = y, perspective = FALSE,
      errors = "bootstrap", bootstrap = "wild", B = 9L, band = "all",
      neval = 6L, legend = 1
    )),
    "legend must be TRUE/FALSE"
  )
})
