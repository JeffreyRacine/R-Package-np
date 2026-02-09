test_that("npconmode basic functionality works", {
  skip_on_cran()
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")

  data("faithful")
  # Use small subset
  # npconmode needs a conditional density bandwidth object
  # where ydat is discrete.
  y <- factor(ifelse(faithful$eruptions[1:50] > 3, "long", "short"))
  x <- faithful$waiting[1:50]
  
  mydat <- data.frame(y, x)
  mpi.bcast.Robj2slave(mydat)

  # bws[1] is for ydat (discrete), bws[2] is for xdat (continuous)
  mpi.bcast.cmd(bw <- npcdensbw(xdat=x, ydat=y, 
                                bws=c(0.1, 5), bandwidth.compute=FALSE),
                caller.execute=TRUE)
  
  mpi.bcast.cmd(mode_fit <- npconmode(bws=bw), caller.execute=TRUE)
  
  expect_s3_class(mode_fit, "conmode")
  expect_type(fitted(mode_fit), "double")
})

test_that("npquantile basic functionality works", {
  skip_on_cran()
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")

  data("faithful")
  x <- faithful$waiting[1:50]
  mpi.bcast.Robj2slave(x)

  mpi.bcast.cmd(bw <- npudistbw(dat=x, bws=5, bandwidth.compute=FALSE),
                caller.execute=TRUE)
  
  # npquantile is for unconditional quantile
  mpi.bcast.cmd(q <- npquantile(x, bws=bw, tau=c(0.25, 0.5, 0.75)),
                caller.execute=TRUE)
  
  expect_type(q, "double")
  expect_length(q, 3)
})
