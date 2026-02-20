test_that("mixed auto->manual plot workflow fails fast", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")

  old.auto <- getOption("npRmpi.autodispatch", FALSE)
  on.exit(options(npRmpi.autodispatch = old.auto), add = TRUE)

  set.seed(41)
  n <- 60
  d <- data.frame(x = rnorm(n), y = rnorm(n))
  mpi.bcast.Robj2slave(d)

  options(npRmpi.autodispatch = TRUE)
  bw <- npregbw(y ~ x, data = d, regtype = "lc", bwmethod = "cv.ls", nmulti = 1)
  fit <- npreg(bws = bw, data = d)

  options(npRmpi.autodispatch = FALSE)
  mpi.bcast.Robj2slave(fit)

  expect_error(
    mpi.bcast.cmd(plot(fit, persp = FALSE, view = "fixed", plot.behavior = "data"),
                  caller.execute = TRUE),
    "cannot be executed inside mpi\\.bcast\\.cmd"
  )
})
