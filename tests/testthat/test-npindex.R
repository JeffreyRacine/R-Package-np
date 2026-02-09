test_that("npindex basic functionality works", {
  skip_on_cran()
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")

  set.seed(42)
  n <- 100
  x1 <- runif(n)
  x2 <- runif(n)
  # Single index model: y = g(x1 + x2) + e
  y <- (x1 + x2)^2 + rnorm(n, sd=0.1)
  
  mydat <- data.frame(y, x1, x2)
  mpi.bcast.Robj2slave(mydat)

  # Optimization might be slow, but for n=100 it should be okay.
  # Using a formula interface
  mpi.bcast.cmd(bw <- npindexbw(y~x1+x2, data=mydat, method="ichimura", nmulti=1),
                caller.execute=TRUE)
  
  mpi.bcast.cmd(model <- npindex(bws=bw), caller.execute=TRUE)
  
  expect_s3_class(model, "singleindex")
  expect_type(predict(model), "double")
  expect_output(summary(model))
})