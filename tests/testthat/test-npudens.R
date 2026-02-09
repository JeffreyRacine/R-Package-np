test_that("npudens basic functionality works", {
  # skip_on_cran()
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")

  data("faithful")
  mpi.bcast.Robj2slave(faithful)

  # Use fixed bandwidths to avoid time-consuming cross-validation in tests
  mpi.bcast.cmd(bw <- npudensbw(dat=faithful, bws=c(0.5, 5), bandwidth.compute=FALSE),
                caller.execute=TRUE)
  
  mpi.bcast.cmd(fit <- npudens(bws=bw), caller.execute=TRUE)
  
  expect_s3_class(fit, "npdensity")
  expect_type(predict(fit), "double")
  expect_equal(length(predict(fit)), nrow(faithful))
  
  # Check summary and print
  expect_output(print(fit))
  expect_output(summary(fit))
})

test_that("npudens works with formula and factors", {
  # skip_on_cran()
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")

  data("Italy")
  Italy_sub <- Italy[1:100, ]
  mpi.bcast.Robj2slave(Italy_sub)

  mpi.bcast.cmd(bw <- npudensbw(formula=~ordered(year)+gdp, data=Italy_sub, 
                                bws=c(0.5, 1.0), bandwidth.compute=FALSE),
                caller.execute=TRUE)
  
  mpi.bcast.cmd(fit <- npudens(bws=bw), caller.execute=TRUE)
  expect_s3_class(fit, "npdensity")
})