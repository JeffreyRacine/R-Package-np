test_that("npksum numeric and formula interfaces agree", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")

  set.seed(1)
  n <- 30
  x <- runif(n)
  y <- rnorm(n)
  dat <- data.frame(y = y, x = x)
  {
    k1 <- npksum(txdat = dat$x, tydat = dat$y, bws = 0.5)
    k2 <- npksum(y ~ x, data = dat, bws = 0.5)
  }

  expect_s3_class(k1, "npkernelsum")
  expect_s3_class(k2, "npkernelsum")
  expect_identical(k1$ksum, k2$ksum)
})

test_that("npksum preserves 1D exdat column naming", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")

  set.seed(1)
  n <- 10
  x <- runif(n)
  ex <- seq(0, 1, length.out = 3)

  {
    k <- npksum(txdat = x, exdat = ex, bws = 0.5)
  }

  expect_identical(colnames(k$eval), "exdat")
})

