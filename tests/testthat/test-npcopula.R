test_that("npcopula basic functionality works", {
  skip_on_cran()
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")

  data("faithful")
  mpi.bcast.Robj2slave(faithful)

  # npcopula can take a dbandwidth (from npudistbw) for copula 
  # or bandwidth (from npudensbw) for copula density
  
  mpi.bcast.cmd(bw <- npudistbw(dat=faithful, bws=c(0.5, 5), bandwidth.compute=FALSE),
                caller.execute=TRUE)
  
  # Copula CDF
  mpi.bcast.cmd(cop <- npcopula(data=faithful, bws=bw), caller.execute=TRUE)
  expect_s3_class(cop, "data.frame")
  expect_true("copula" %in% names(cop))
  expect_true(all(cop$copula >= 0 & cop$copula <= 1))
  
  # Copula density
  mpi.bcast.cmd(bw_dens <- npudensbw(dat=faithful, bws=c(0.5, 5), bandwidth.compute=FALSE),
                caller.execute=TRUE)
  mpi.bcast.cmd(cop_dens <- npcopula(data=faithful, bws=bw_dens), caller.execute=TRUE)
  expect_s3_class(cop_dens, "data.frame")
  expect_true("copula" %in% names(cop_dens))
  expect_true(all(cop_dens$copula >= 0))
})