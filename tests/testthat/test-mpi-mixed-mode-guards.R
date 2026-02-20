test_that("autodispatch plot workflow stays consistent", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")

  old.auto <- getOption("npRmpi.autodispatch", FALSE)
  on.exit(options(npRmpi.autodispatch = old.auto), add = TRUE)

  set.seed(41)
  n <- 60
  d <- data.frame(x = rnorm(n), y = rnorm(n))

  options(npRmpi.autodispatch = TRUE)
  bw <- npregbw(y ~ x, data = d, regtype = "lc", bwmethod = "cv.ls", nmulti = 1)
  fit <- npreg(bws = bw, data = d)

  out <- plot(fit, persp = FALSE, view = "fixed", plot.behavior = "data")
  expect_true(is.list(out))
})
