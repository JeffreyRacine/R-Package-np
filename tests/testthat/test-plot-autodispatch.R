test_that("plot works with autodispatch for non-bootstrap paths", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")

  old.auto <- getOption("npRmpi.autodispatch", FALSE)
  on.exit(options(npRmpi.autodispatch = old.auto), add = TRUE)
  options(npRmpi.autodispatch = TRUE)

  set.seed(21)
  n <- 60
  d <- data.frame(x = rnorm(n), y = rnorm(n))

  bw <- npregbw(y ~ x, data = d, regtype = "lc", bwmethod = "cv.ls", nmulti = 1)
  out <- suppressWarnings(plot(bw,
                               persp = FALSE,
                               view = "fixed",
                               plot.behavior = "data",
                               plot.errors.method = "none"))

  expect_type(out, "list")
  expect_true(length(out) > 0)
})

test_that("plot bootstrap path works under autodispatch", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")

  old.auto <- getOption("npRmpi.autodispatch", FALSE)
  on.exit(options(npRmpi.autodispatch = old.auto), add = TRUE)
  options(npRmpi.autodispatch = TRUE)

  set.seed(22)
  n <- 60
  d <- data.frame(x = rnorm(n), y = rnorm(n))

  bw <- npregbw(y ~ x, data = d, regtype = "lc", bwmethod = "cv.ls", nmulti = 1)

  out <- suppressWarnings(
    plot(bw,
         persp = FALSE,
         view = "fixed",
         plot.behavior = "data",
         plot.errors.method = "bootstrap",
         plot.errors.boot.num = 9)
  )

  expect_type(out, "list")
  expect_true(length(out) > 0)
})
