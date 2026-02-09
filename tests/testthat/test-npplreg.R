test_that("npplreg basic functionality works", {
  skip_on_cran()
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")

  set.seed(42)
  n <- 100
  x1 <- runif(n) # nonparametric part
  z1 <- runif(n) # parametric part
  y <- x1^2 + 2*z1 + rnorm(n, sd=0.1)
  
  mydat <- data.frame(y, x1, z1)
  mpi.bcast.Robj2slave(mydat)

  # Partially linear model: y = g(x1) + z1*beta + e
  # bws needs to be a matrix. Row 1: y on x1, Row 2: z1 on x1
  mpi.bcast.cmd(bw_mat <- matrix(c(0.1, 0.1), nrow=2, ncol=1), caller.execute=TRUE)
  mpi.bcast.cmd(bw <- npplregbw(xdat=z1, zdat=x1, ydat=y, bws=bw_mat, bandwidth.compute=FALSE),
                caller.execute=TRUE)
  
  mpi.bcast.cmd(model <- npplreg(bws=bw), caller.execute=TRUE)
  
  expect_s3_class(model, "plregression")
  expect_type(predict(model), "double")
  expect_output(summary(model))
  expect_length(coef(model), 1)
})