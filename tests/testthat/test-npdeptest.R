test_that("npdeptest basic functionality works", {
  skip_on_cran()
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")

  set.seed(42)
  n <- 50
  x <- rnorm(n)
  y <- x + rnorm(n, sd=0.1)
  
  mpi.bcast.Robj2slave(x)
  mpi.bcast.Robj2slave(y)

  # Use summation and small boot.num for speed
  mpi.bcast.cmd(test <- npdeptest(x, y, method="summation", boot.num=19),
                caller.execute=TRUE)
  
  expect_s3_class(test, "deptest")
  expect_output(summary(test))
})

test_that("npsdeptest basic functionality works", {
  skip_on_cran()
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")

  set.seed(42)
  n <- 50
  y <- arima.sim(n=n, list(ar=0.5))
  mpi.bcast.Robj2slave(y)
  
  mpi.bcast.cmd(test <- npsdeptest(y, lag.num=1, method="summation", boot.num=19),
                caller.execute=TRUE)
  
  expect_s3_class(test, "sdeptest")
  expect_output(summary(test))
})