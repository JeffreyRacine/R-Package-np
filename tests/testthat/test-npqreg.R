test_that("npqreg basic functionality works", {
  # skip_on_cran()
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")

  data("cps71")
  cps71_sub <- cps71[1:50, ]
  mpi.bcast.Robj2slave(cps71_sub)
  
  # Quantile regression needs a condbandwidth object
  mpi.bcast.cmd(bw <- npcdistbw(xdat=cps71_sub$age, ydat=cps71_sub$logwage, 
                                bws=c(0.5, 5), bandwidth.compute=FALSE),
                caller.execute=TRUE)
  
  # Median regression
  mpi.bcast.cmd(model <- npqreg(bws=bw, tau=0.5), caller.execute=TRUE)
  
  expect_s3_class(model, "qregression")
  expect_type(predict(model), "double")
  expect_equal(length(predict(model)), 50)
  expect_equal(model$tau, 0.5)
  
  expect_output(summary(model))
})

test_that("npqreg works with multiple taus", {
  # skip_on_cran()
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")

  data("cps71")
  cps71_sub <- cps71[1:30, ]
  mpi.bcast.Robj2slave(cps71_sub)

  mpi.bcast.cmd(bw <- npcdistbw(xdat=cps71_sub$age, ydat=cps71_sub$logwage, 
                                bws=c(0.5, 5), bandwidth.compute=FALSE),
                caller.execute=TRUE)
  
  mpi.bcast.cmd(model_q25 <- npqreg(bws=bw, tau=0.25), caller.execute=TRUE)
  expect_equal(model_q25$tau, 0.25)
})