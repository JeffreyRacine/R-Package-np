test_that("npudist basic functionality works", {
  # skip_on_cran()
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")

  data("faithful")
  mpi.bcast.Robj2slave(faithful)

  # Use fixed bandwidths for speed
  mpi.bcast.cmd(bw <- npudistbw(dat=faithful, bws=c(0.5, 5), bandwidth.compute=FALSE),
                caller.execute=TRUE)
  
  mpi.bcast.cmd(fit <- npudist(bws=bw), caller.execute=TRUE)
  
  expect_s3_class(fit, "npdistribution")
  expect_type(predict(fit), "double")
  expect_true(all(predict(fit) >= 0 & predict(fit) <= 1))
  expect_equal(length(predict(fit)), nrow(faithful))
  
  expect_output(summary(fit))
})

test_that("npudist works with formula and factors", {
  # skip_on_cran()
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")

  data("Italy")
  Italy_sub <- Italy[1:100, ]
  mpi.bcast.Robj2slave(Italy_sub)

  mpi.bcast.cmd(bw <- npudistbw(formula=~ordered(year)+gdp, data=Italy_sub, 
                                bws=c(0.5, 1.0), bandwidth.compute=FALSE),
                caller.execute=TRUE)
  
  mpi.bcast.cmd(fit <- npudist(bws=bw), caller.execute=TRUE)
  expect_s3_class(fit, "npdistribution")
})