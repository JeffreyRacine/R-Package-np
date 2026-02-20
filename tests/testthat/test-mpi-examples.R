context("MPI Examples")

test_that("MPI examples from man pages work correctly", {
  # Skip if Rmpi is not functional in this environment
  # skip_on_cran()
  
  if (!spawn_mpi_slaves()) {
    skip("Could not spawn MPI slaves for testing")
  }

  # Ensure cleanup on exit
  on.exit(try(close_mpi_slaves(force=TRUE), silent=TRUE))
  options(npRmpi.autodispatch = TRUE, np.messages = FALSE)

  # 1. npregbw example (from npseed.Rd correction)
  test_that("npregbw works in parallel", {
    npseed(712)
    
    n <- 50 # Small n for fast test
    x <- runif(n)
    y <- x + rnorm(n, sd = 0.1)
    mydat <- data.frame(x,y)
    
    bw <- npregbw(y~x, data=mydat)
    
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
    
    model.iv <- with(Engel_sub, npregiv(y=food, z=logexp, w=logwages,
                                        method="Landweber-Fridman", iterate.max=2))
    
    expect_s3_class(model.iv, "npregiv")
  })

  # 3. np.pairs example
  test_that("np.pairs works in parallel", {
    data("USArrests")
    # Small subset
    dat <- USArrests[1:20, ]
    y_vars <- c("Murder", "UrbanPop")
    
    pair_list <- np.pairs(y_vars = y_vars, y_dat = dat,
                          bandwidth.compute=FALSE, bws=0.5)
    
    expect_type(pair_list, "list")
    expect_equal(length(pair_list$pair_kerns), 4)
  })

})
