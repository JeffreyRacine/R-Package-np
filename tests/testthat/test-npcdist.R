test_that("npcdist basic functionality works", {
  skip_on_cran()
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")

  data("faithful")
  mpi.bcast.Robj2slave(faithful)

  mpi.bcast.cmd(bw <- npcdistbw(xdat=faithful$eruptions, ydat=faithful$waiting, 
                                bws=c(0.5, 5), bandwidth.compute=FALSE),
                caller.execute=TRUE)
  
  mpi.bcast.cmd(fit <- npcdist(bws=bw), caller.execute=TRUE)
  
  expect_s3_class(fit, "condistribution")
  expect_type(predict(fit), "double")
  expect_equal(length(predict(fit)), nrow(faithful))
  
  expect_output(summary(fit))
})