test_that("npregiv basic functionality works", {
  # skip_on_cran()
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")

  set.seed(42)
  n <- 100
  # Simple IV setup
  z <- runif(n) # instrument
  w <- z + rnorm(n, sd=0.1) # endogenous regressor
  y <- w^2 + rnorm(n, sd=0.1) # outcome
  
  mydat <- data.frame(y, w, z)
  mpi.bcast.Robj2slave(mydat)

  # Let's try a very small example
  mpi.bcast.cmd(model <- npregiv(y=y, w=w, z=z, method="Landweber-Fridman", iterate.max=10),
                caller.execute=TRUE)
  
  expect_s3_class(model, "npregiv")
  expect_true("phi" %in% names(model))
  expect_false(is.null(model$call))
})

test_that("npregivderiv basic functionality works", {
  # skip_on_cran()
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")

  set.seed(42)
  n <- 50
  z <- runif(n)
  w <- z + rnorm(n, sd=0.1)
  y <- z^2 + rnorm(n, sd=0.1)
  
  mydat <- data.frame(y, z, w)
  mpi.bcast.Robj2slave(mydat)

  # Just check if it runs
  mpi.bcast.cmd(model <- npregivderiv(y=y, z=z, w=w, iterate.max=2),
                caller.execute=TRUE)
  expect_s3_class(model, "npregivderiv")
  expect_true("phi.prime" %in% names(model))
  expect_false(is.null(model$call))
})