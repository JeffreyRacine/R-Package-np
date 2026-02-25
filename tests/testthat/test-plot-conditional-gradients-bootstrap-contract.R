test_that("conditional density/distribution gradient bootstrap inid works for bandwidth objects", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  old.auto <- getOption("npRmpi.autodispatch", FALSE)
  on.exit(options(npRmpi.autodispatch = old.auto), add = TRUE)
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  options(npRmpi.autodispatch = TRUE)

  set.seed(20260225)
  n <- 40L
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  y <- rnorm(n)

  xdat <- data.frame(x1 = x1, x2 = x2)
  ydat <- data.frame(y = y)

  bw.cd <- npcdensbw(
    xdat = xdat,
    ydat = ydat,
    bws = c(0.45, 0.45, 0.45),
    bandwidth.compute = FALSE
  )
  bw.cdist <- npcdistbw(
    xdat = xdat,
    ydat = ydat,
    bws = c(0.45, 0.45, 0.45),
    bandwidth.compute = FALSE
  )

  out.cd <- suppressWarnings(
    plot(
      bw.cd,
      plot.behavior = "data",
      perspective = FALSE,
      gradients = TRUE,
      plot.errors.method = "bootstrap",
      plot.errors.boot.method = "inid",
      plot.errors.boot.num = 5
    )
  )
  out.cdist <- suppressWarnings(
    plot(
      bw.cdist,
      plot.behavior = "data",
      perspective = FALSE,
      gradients = TRUE,
      plot.errors.method = "bootstrap",
      plot.errors.boot.method = "inid",
      plot.errors.boot.num = 5
    )
  )

  expect_type(out.cd, "list")
  expect_type(out.cdist, "list")
  expect_true(length(out.cd) >= 1L)
  expect_true(length(out.cdist) >= 1L)
  expect_true(all(vapply(out.cd, function(z) "gc1err" %in% names(z), logical(1))))
  expect_true(all(vapply(out.cdist, function(z) "gc1err" %in% names(z), logical(1))))
})

test_that("conditional density/distribution gradient bootstrap inid works for fitted objects", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  old.auto <- getOption("npRmpi.autodispatch", FALSE)
  on.exit(options(npRmpi.autodispatch = old.auto), add = TRUE)
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  options(npRmpi.autodispatch = TRUE)

  set.seed(20260226)
  n <- 35L
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  y <- rnorm(n)

  xdat <- data.frame(x1 = x1, x2 = x2)
  ydat <- data.frame(y = y)

  bw.cd <- npcdensbw(
    xdat = xdat,
    ydat = ydat,
    bws = c(0.5, 0.5, 0.5),
    bandwidth.compute = FALSE
  )
  bw.cdist <- npcdistbw(
    xdat = xdat,
    ydat = ydat,
    bws = c(0.5, 0.5, 0.5),
    bandwidth.compute = FALSE
  )

  fit.cd <- npcdens(bws = bw.cd)
  fit.cdist <- npcdist(bws = bw.cdist)

  out.cd <- suppressWarnings(
    plot(
      fit.cd,
      plot.behavior = "data",
      perspective = FALSE,
      gradients = TRUE,
      plot.errors.method = "bootstrap",
      plot.errors.boot.method = "inid",
      plot.errors.boot.num = 5
    )
  )
  out.cdist <- suppressWarnings(
    plot(
      fit.cdist,
      plot.behavior = "data",
      perspective = FALSE,
      gradients = TRUE,
      plot.errors.method = "bootstrap",
      plot.errors.boot.method = "inid",
      plot.errors.boot.num = 5
    )
  )

  expect_type(out.cd, "list")
  expect_type(out.cdist, "list")
  expect_true(length(out.cd) >= 1L)
  expect_true(length(out.cdist) >= 1L)
})
