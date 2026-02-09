context("MPI Examples")

test_that("MPI examples from man pages work correctly", {
  # Skip if Rmpi is not functional in this environment
  # skip_on_cran()
  
  if (!requireNamespace("Rmpi", quietly = TRUE)) {
    skip("Rmpi package not available")
  }

  # Attempt to spawn slaves
  # Note: In some test environments, spawning might fail depending on MPI config
  # We wrap in try to be safe
  spawn_status <- try(mpi.spawn.Rslaves(nslaves=1, quiet=TRUE), silent=TRUE)
  
  if (inherits(spawn_status, "try-error") || mpi.comm.size(0) < 2) {
    skip("Could not spawn MPI slaves for testing")
  }

  # Ensure cleanup on exit
  on.exit(try(mpi.close.Rslaves(), silent=TRUE))

  # Initialize
  mpi.bcast.cmd(np.mpi.initialize(), caller.execute=TRUE)

  # 1. npregbw example (from npseed.Rd correction)
  test_that("npregbw works in parallel", {
    mpi.bcast.cmd(npseed(712), caller.execute=TRUE)
    
    n <- 50 # Small n for fast test
    x <- runif(n)
    y <- x + rnorm(n, sd = 0.1)
    mydat <- data.frame(x,y)
    
    mpi.bcast.Robj2slave(mydat)
    
    mpi.bcast.cmd(bw <- npregbw(y~x, data=mydat), caller.execute=TRUE)
    
    # We can retrieve the result if needed, or just check it exists on slaves
    # But usually we check it on master
    expect_s3_class(bw, "rbandwidth")
  })

  # 2. npregiv example (from Engel95.Rd)
  test_that("npregiv works in parallel", {
    data(Engel95)
    # Use small subset for testing
    Engel_sub <- Engel95[1:100, ]
    Engel_sub <- Engel_sub[order(Engel_sub$logexp),]
    
    mpi.bcast.Robj2slave(Engel_sub)
    mpi.bcast.cmd(attach(Engel_sub), caller.execute=TRUE)
    
    mpi.bcast.cmd(model.iv <- npregiv(y=food, z=logexp, w=logwages, 
                                     method="Landweber-Fridman", iterate.max=2),
                  caller.execute=TRUE)
    
    expect_s3_class(model.iv, "npregiv")
    
    mpi.bcast.cmd(detach(Engel_sub), caller.execute=TRUE)
  })

  # 3. np.pairs example
  test_that("np.pairs works in parallel", {
    data("USArrests")
    # Small subset
    dat <- USArrests[1:20, ]
    y_vars <- c("Murder", "UrbanPop")
    
    mpi.bcast.Robj2slave(dat)
    mpi.bcast.Robj2slave(y_vars)
    
    mpi.bcast.cmd(pair_list <- np.pairs(y_vars = y_vars, y_dat = dat, 
                                        bandwidth.compute=FALSE, bws=0.5),
                  caller.execute=TRUE)
    
    expect_type(pair_list, "list")
    expect_equal(length(pair_list$pair_kerns), 4)
  })

})
