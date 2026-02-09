test_that("npscoef basic functionality works", {
  # skip_on_cran()
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")

  set.seed(42)
  n <- 100
  x1 <- runif(n) # smoothing variable
  z1 <- runif(n) # parametric variable
  # Smooth coefficient model: y = beta(x1) * z1 + e
  y <- (x1^2) * z1 + rnorm(n, sd=0.1)
  
  mydat <- data.frame(y, x1, z1)
  mpi.bcast.Robj2slave(mydat)

  mpi.bcast.cmd(bw <- npscoefbw(xdat=x1, zdat=z1, ydat=y, bws=0.1, bandwidth.compute=FALSE),
                caller.execute=TRUE)
  
  mpi.bcast.cmd(model <- npscoef(bws=bw), caller.execute=TRUE)
  
  expect_s3_class(model, "smoothcoefficient")
  expect_type(predict(model), "double")
  expect_output(summary(model))
})