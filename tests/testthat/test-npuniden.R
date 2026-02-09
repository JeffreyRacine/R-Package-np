test_that("npuniden.boundary basic functionality works", {
  skip_on_cran()
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")

  set.seed(42)
  n <- 50
  X <- rbeta(n, 5, 1)
  mpi.bcast.Robj2slave(X)
  
  # Use fixed bandwidth for speed
  mpi.bcast.cmd(model <- npuniden.boundary(X, h=0.05, a=0, b=1), caller.execute=TRUE)
  
  expect_type(model, "list")
  expect_true("f" %in% names(model))
  expect_length(model$f, n)
  expect_true(all(model$f >= 0))
})

test_that("npuniden.reflect basic functionality works", {
  skip_on_cran()
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")

  set.seed(42)
  n <- 50
  X <- rbeta(n, 5, 1)
  mpi.bcast.Robj2slave(X)
  
  mpi.bcast.cmd(model <- npuniden.reflect(X, h=0.05, a=0, b=1), caller.execute=TRUE)
  
  expect_type(model, "list")
  expect_true("f" %in% names(model))
})

test_that("npuniden.sc basic functionality works", {
  skip_on_cran()
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")

  set.seed(42)
  n <- 50
  X <- rbeta(n, 5, 1)
  mpi.bcast.Robj2slave(X)
  
  # Shape constrained density estimation
  mpi.bcast.cmd(model <- npuniden.sc(X, h=0.05, a=0, b=1, lb=0, ub=Inf, constraint="density"),
                caller.execute=TRUE)
  
  expect_type(model, "list")
  expect_true("f" %in% names(model))
})