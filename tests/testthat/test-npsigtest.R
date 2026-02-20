test_that("npsigtest basic functionality works", {
  # skip_on_cran()
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  old.auto <- getOption("npRmpi.autodispatch", FALSE)
  options(npRmpi.autodispatch = FALSE)
  on.exit(options(npRmpi.autodispatch = old.auto), add = TRUE)
  mpi.bcast.cmd(options(npRmpi.autodispatch = FALSE), caller.execute = TRUE)

  set.seed(42)
  n <- 50 # Keep it small for speed
  x1 <- runif(n)
  x2 <- runif(n)
  y <- x1^2 + rnorm(n, sd=0.1) # x2 is irrelevant
  
  mydat <- data.frame(y, x1, x2)
  mpi.bcast.Robj2slave(mydat)

  mpi.bcast.cmd(bw <- npregbw(y~x1+x2, data=mydat, bws=c(0.1, 0.5), bandwidth.compute=FALSE),
                caller.execute=TRUE)
  
  # Significance test can be slow, use few boot replications
  mpi.bcast.cmd(sig <- npsigtest(bws=bw, boot.num=19),
                caller.execute=TRUE)
  
  expect_s3_class(sig, "sigtest")
  expect_output(summary(sig))
})
