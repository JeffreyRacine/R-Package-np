library(npRmpi)

test_that("plot condensity 2D data payload applies proper repair on supported grids", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")

  set.seed(2)
  n <- 80L
  x <- runif(n, -1, 1)
  y <- sin(2 * pi * x) + rnorm(n, sd = 0.20)

  bw <- npcdensbw(
    xdat = data.frame(x = x),
    ydat = data.frame(y = y),
    bws = c(0.28, 0.22),
    bandwidth.compute = FALSE,
    regtype = "lp",
    degree = 3L
  )
  fit <- npcdens(
    bws = bw,
    txdat = data.frame(x = x),
    tydat = data.frame(y = y),
    proper = TRUE
  )

  out <- plot(
    fit,
    plot.behavior = "data",
    perspective = TRUE,
    view = "fixed"
  )

  expect_type(out, "list")
  expect_true("cd1" %in% names(out))
  expect_s3_class(out$cd1, "condensity")
  expect_true(isTRUE(out$cd1$proper.requested))
  expect_true(isTRUE(out$cd1$proper.applied))
  expect_true(all(out$cd1$condens >= -1e-8))
})

test_that("plot condensity 1D data payload repairs only y-varying panels", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")

  set.seed(2)
  n <- 80L
  x <- runif(n, -1, 1)
  y <- sin(2 * pi * x) + rnorm(n, sd = 0.20)

  bw <- npcdensbw(
    xdat = data.frame(x = x),
    ydat = data.frame(y = y),
    bws = c(0.28, 0.22),
    bandwidth.compute = FALSE,
    regtype = "lp",
    degree = 3L
  )
  fit <- npcdens(
    bws = bw,
    txdat = data.frame(x = x),
    tydat = data.frame(y = y),
    proper = TRUE
  )

  out <- suppressWarnings(plot(
    fit,
    plot.behavior = "data",
    perspective = FALSE
  ))

  expect_length(out, 2L)
  expect_s3_class(out[[1]], "condensity")
  expect_s3_class(out[[2]], "condensity")
  expect_true(isTRUE(out[[1]]$proper.requested))
  expect_false(isTRUE(out[[1]]$proper.applied))
  expect_true(isTRUE(out[[2]]$proper.requested))
  expect_true(isTRUE(out[[2]]$proper.applied))
  expect_true(all(out[[2]]$condens >= -1e-8))
})

test_that("plot condensity rejects asymptotic errors when repair is active", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")

  set.seed(2)
  x <- runif(60, -1, 1)
  y <- sin(2 * pi * x) + rnorm(60, sd = 0.2)

  bw <- npcdensbw(
    xdat = data.frame(x = x),
    ydat = data.frame(y = y),
    bws = c(0.28, 0.22),
    bandwidth.compute = FALSE,
    regtype = "lp",
    degree = 3L
  )
  fit <- npcdens(
    bws = bw,
    txdat = data.frame(x = x),
    tydat = data.frame(y = y),
    proper = TRUE
  )

  expect_error(
    suppressWarnings(plot(
      fit,
      plot.behavior = "data",
      perspective = TRUE,
      view = "fixed",
      plot.errors.method = "asymptotic"
    )),
    "unsupported when proper=TRUE"
  )
})
