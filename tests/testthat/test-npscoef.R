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
  bw <- npscoefbw(xdat=x1, zdat=z1, ydat=y, bws=0.1, bandwidth.compute=FALSE)
  
  model <- npscoef(bws=bw)
  
  expect_s3_class(model, "smoothcoefficient")
  expect_type(predict(model), "double")
  expect_output(summary(model))
})

test_that("npscoefbw records ll/lp controls", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  set.seed(43)
  n <- 80
  x1 <- runif(n)
  z1 <- runif(n)
  y <- (0.2 + x1) * cos(2 * pi * z1) + rnorm(n, sd = 0.08)

  bw.ll <- npscoefbw(
    xdat = x1,
    zdat = z1,
    ydat = y,
    regtype = "ll",
    bws = 0.15,
    bandwidth.compute = FALSE
  )
  expect_identical(bw.ll$regtype, "ll")

  bw.lp <- npscoefbw(
    xdat = x1,
    zdat = z1,
    ydat = y,
    regtype = "lp",
    basis = "tensor",
    degree = 2L,
    bws = 0.15,
    bandwidth.compute = FALSE
  )
  expect_identical(bw.lp$regtype, "lp")
  expect_identical(bw.lp$basis, "tensor")
})
