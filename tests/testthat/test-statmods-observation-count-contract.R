suppressPackageStartupMessages(library(npRmpi))

test_that("npRmpi nonfixed search accepts observation-count k beyond empirical support size", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")

  x <- data.frame(x = c(0, 0, 0, 1))
  y <- c(0, 1, 2, 3)

  bw.adap <- npregbw(
    ydat = y,
    xdat = x,
    bwtype = "adaptive_nn",
    bws = length(y) - 1L,
    bandwidth.compute = FALSE
  )
  bw.gen <- npregbw(
    ydat = y,
    xdat = x,
    bwtype = "generalized_nn",
    bws = 2L,
    bandwidth.compute = FALSE
  )

  expect_no_error(
    npreg(bws = bw.adap, exdat = x)
  )
  expect_no_error(
    npreg(bws = bw.gen, exdat = x)
  )
})

test_that("npRmpi generalized observation-count radii handle low-support ties", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")

  x <- data.frame(x = c(0, 0, 0, 1, 1, 2))
  y <- c(0, 1, 2, 3, 4, 5)
  ex <- data.frame(x = c(0, 0.5, 1, 1.5, 2))

  bw <- npregbw(
    ydat = y,
    xdat = x,
    bwtype = "generalized_nn",
    bws = 3L,
    bandwidth.compute = FALSE
  )

  fit <- npreg(bws = bw, exdat = ex)

  expect_length(fit$mean, nrow(ex))
  expect_true(all(is.finite(fit$mean)))
})
