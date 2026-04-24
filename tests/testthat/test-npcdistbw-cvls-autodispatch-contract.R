library(npRmpi)

test_that("npcdistbw cv.ls autodispatch does not require npcdens quadrature slots", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  set.seed(20260424)
  n <- 30L
  x <- rnorm(n)
  y <- 0.4 * x + rnorm(n)
  dat <- data.frame(x = x, y = y)

  bw <- npcdistbw(
    y ~ x,
    data = dat,
    bwmethod = "cv.ls",
    nmulti = 1L,
    itmax = 1L
  )

  expect_identical(bw$method, "cv.ls")
  expect_true(is.finite(bw$fval))
})
