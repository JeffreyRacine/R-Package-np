test_that("npcdens basic functionality works", {
  skip_on_cran()
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")

  data("faithful")
  mpi.bcast.Robj2slave(faithful)

  # waiting is Y, eruptions is X
  mpi.bcast.cmd(bw <- npcdensbw(xdat=faithful$eruptions, ydat=faithful$waiting, 
                                bws=c(0.5, 5), bandwidth.compute=FALSE),
                caller.execute=TRUE)
  
  mpi.bcast.cmd(fit <- npcdens(bws=bw), caller.execute=TRUE)
  
  expect_s3_class(fit, "condensity")
  expect_type(predict(fit), "double")
  expect_equal(length(predict(fit)), nrow(faithful))
  
  expect_output(summary(fit))
})

test_that("npcdens works with formula", {
  skip_on_cran()
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")

  data("Italy")
  Italy_sub <- Italy[1:50, ]
  mpi.bcast.Robj2slave(Italy_sub)

  mpi.bcast.cmd(bw <- npcdensbw(gdp~ordered(year), data=Italy_sub, 
                                bws=c(0.5, 1.0), bandwidth.compute=FALSE),
                caller.execute=TRUE)
  
  mpi.bcast.cmd(fit <- npcdens(bws=bw), caller.execute=TRUE)
  expect_s3_class(fit, "condensity")
})