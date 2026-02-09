test_that("npcmstest basic functionality works", {
  skip_on_cran()
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")

  set.seed(42)
  n <- 50
  x <- rnorm(n)
  y <- 1 + x + rnorm(n, sd=0.1)
  
  mydat <- data.frame(x, y)
  mpi.bcast.Robj2slave(mydat)

  # Parametric model - need x=TRUE, y=TRUE for npcmstest
  mpi.bcast.cmd(model <- lm(y~x, data=mydat, x=TRUE, y=TRUE), caller.execute=TRUE)
  
  # Asymptotic for speed
  mpi.bcast.cmd(test <- npcmstest(model=model, xdat=x, ydat=y, distribution="asymptotic"),
                caller.execute=TRUE)
  
  expect_s3_class(test, "cmstest")
  expect_output(summary(test))
})

test_that("npqcmstest basic functionality works", {
  skip_on_cran()
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")

  mpi.bcast.cmd(library(quantreg), caller.execute=TRUE)
  set.seed(42)
  n <- 50
  x <- rnorm(n)
  y <- 1 + x + rnorm(n, sd=0.1)
  
  mydat <- data.frame(x, y)
  mpi.bcast.Robj2slave(mydat)

  # npqcmstest needs an rq model with model=TRUE
  mpi.bcast.cmd(model <- rq(y~x, data=mydat, tau=0.5, model=TRUE), caller.execute=TRUE)
  
  # Asymptotic for speed
  mpi.bcast.cmd(test <- npqcmstest(model=model, xdat=x, ydat=y, distribution="asymptotic"),
                caller.execute=TRUE)
  
  expect_s3_class(test, "cmstest") 
})

test_that("npdeneqtest basic functionality works", {
  skip_on_cran()
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")

  set.seed(42)
  n <- 50
  x <- data.frame(v1=rnorm(n))
  y <- data.frame(v1=rnorm(n, mean=0.5))
  
  mpi.bcast.Robj2slave(x)
  mpi.bcast.Robj2slave(y)

  # Use small boot.num
  mpi.bcast.cmd(test <- npdeneqtest(x, y, boot.num=19), caller.execute=TRUE)
  
  expect_s3_class(test, "deneqtest")
  expect_output(summary(test))
})

test_that("npsymtest basic functionality works", {
  skip_on_cran()
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")

  set.seed(42)
  n <- 50
  x <- rgamma(n, shape=2)
  mpi.bcast.Robj2slave(x)
  
  mpi.bcast.cmd(test <- npsymtest(x, method="summation", boot.num=19),
                caller.execute=TRUE)
  
  expect_s3_class(test, "symtest")
  expect_output(summary(test))
})

test_that("npunitest basic functionality works", {
  skip_on_cran()
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")

  set.seed(42)
  n <- 50
  x <- rnorm(n)
  y <- rnorm(n, mean=0.5)
  
  mpi.bcast.Robj2slave(x)
  mpi.bcast.Robj2slave(y)

  mpi.bcast.cmd(test <- npunitest(x, y, method="summation", boot.num=19),
                caller.execute=TRUE)
  
  expect_s3_class(test, "unitest")
  expect_output(summary(test))
})