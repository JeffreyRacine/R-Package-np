context("MPI Comprehensive Examples")

test_that("Comprehensive MPI examples work correctly", {
  skip_on_cran()
  
  if (!requireNamespace("Rmpi", quietly = TRUE)) {
    skip("Rmpi package not available")
  }

  # Attempt to spawn slaves
  spawn_status <- try(mpi.spawn.Rslaves(nslaves=1, quiet=TRUE), silent=TRUE)
  
  if (inherits(spawn_status, "try-error") || mpi.comm.size(0) < 2) {
    skip("Could not spawn MPI slaves for testing")
  }

  # Ensure cleanup on exit
  on.exit(try(mpi.close.Rslaves(), silent=TRUE))

  # Initialize
  mpi.bcast.cmd(np.mpi.initialize(), caller.execute=TRUE)

  # 1. npudensbw / npudens (Unconditional density)
  test_that("npudens works in parallel", {
    n <- 50
    x <- rnorm(n)
    mpi.bcast.Robj2slave(x)
    mpi.bcast.cmd(bw <- npudensbw(~x, bwmethod="normal-reference"), caller.execute=TRUE)
    mpi.bcast.cmd(fhat <- npudens(bw), caller.execute=TRUE)
    expect_s3_class(fhat, "npdensity")
  })

  # 2. npcdensbw / npcdens (Conditional density)
  test_that("npcdens works in parallel", {
    data("Italy")
    Italy_sub <- Italy[1:50, ]
    mpi.bcast.Robj2slave(Italy_sub)
    mpi.bcast.cmd(bw <- npcdensbw(gdp~ordered(year), data=Italy_sub, bwmethod="normal-reference"), caller.execute=TRUE)
    mpi.bcast.cmd(fhat <- npcdens(bw), caller.execute=TRUE)
    expect_s3_class(fhat, "condensity")
  })

  # 3. npudistbw / npudist (Unconditional distribution)
  test_that("npudist works in parallel", {
    n <- 50
    x <- rnorm(n)
    mpi.bcast.Robj2slave(x)
    mpi.bcast.cmd(bw <- npudistbw(~x, bwmethod="normal-reference"), caller.execute=TRUE)
    mpi.bcast.cmd(Fhat <- npudist(bw), caller.execute=TRUE)
    expect_s3_class(Fhat, "npdistribution")
  })

  # 4. npcdistbw / npcdist (Conditional distribution)
  test_that("npcdist works in parallel", {
    data("Italy")
    Italy_sub <- Italy[1:50, ]
    mpi.bcast.Robj2slave(Italy_sub)
    mpi.bcast.cmd(bw <- npcdistbw(gdp~ordered(year), data=Italy_sub, bwmethod="normal-reference"), caller.execute=TRUE)
    mpi.bcast.cmd(Fhat <- npcdist(bw), caller.execute=TRUE)
    expect_s3_class(Fhat, "condistribution")
  })

  # 5. npplregbw / npplreg (Partially linear regression)
  test_that("npplreg works in parallel", {
    n <- 50
    x1 <- rnorm(n)
    x2 <- factor(rbinom(n, 1, .5))
    z1 <- factor(rbinom(n, 1, .5))
    z2 <- rnorm(n)
    y <- 1 + x1 + as.numeric(x2) + as.numeric(z1) + sin(z2) + rnorm(n)
    mydat <- data.frame(y,x1,x2,z1,z2)
    mpi.bcast.Robj2slave(mydat)
    mpi.bcast.cmd(bw <- npplregbw(formula=y~x1+x2|z1+z2, data=mydat, bwmethod="normal-reference"), caller.execute=TRUE)
    mpi.bcast.cmd(pl <- npplreg(bws=bw), caller.execute=TRUE)
    expect_s3_class(pl, "plregression")
  })

  # 6. npindexbw / npindex (Single index model)
  test_that("npindex works in parallel", {
    n <- 50
    x1 <- runif(n, min=-1, max=1)
    x2 <- runif(n, min=-1, max=1)
    y <- x1 - x2 + rnorm(n)
    mydat <- data.frame(x1,x2,y)
    mpi.bcast.Robj2slave(mydat)
    # Use small nmulti for speed
    mpi.bcast.cmd(bw <- npindexbw(formula=y~x1+x2, data=mydat, nmulti=1), caller.execute=TRUE)
    mpi.bcast.cmd(model <- npindex(bws=bw), caller.execute=TRUE)
    expect_s3_class(model, "singleindex")
  })

  # 7. npscoefbw / npscoef (Smooth coefficient)
  test_that("npscoef works in parallel", {
    n <- 50
    x <- runif(n)
    z <- runif(n, min=-2, max=2)
    y <- x*exp(z)*(1.0+rnorm(n,sd = 0.2))
    mydat <- data.frame(x,y,z)
    mpi.bcast.Robj2slave(mydat)
    mpi.bcast.cmd(bw <- npscoefbw(y~x|z, data=mydat, bwmethod="normal-reference"), caller.execute=TRUE)
    mpi.bcast.cmd(model <- npscoef(bws=bw), caller.execute=TRUE)
    expect_s3_class(model, "smoothcoefficient")
  })

  # 8. npsigtest (Significance test)
  test_that("npsigtest works in parallel", {
    n <- 50
    z <- factor(rbinom(n,1,.5))
    x1 <- rnorm(n)
    y <- x1 + rnorm(n)
    mydat <- data.frame(z,x1,y)
    mpi.bcast.Robj2slave(mydat)
    mpi.bcast.cmd(model <- npreg(y~z+x1, data=mydat, bwmethod="normal-reference"), caller.execute=TRUE)
    mpi.bcast.cmd(output <- npsigtest(model, boot.num=2), caller.execute=TRUE)
    expect_s3_class(output, "sigtest")
  })

  # 9. npdeptest / npsdeptest
  test_that("npdeptest works in parallel", {
    n <- 50
    x <- rnorm(n)
    y <- x + rnorm(n)
    mpi.bcast.Robj2slave(x)
    mpi.bcast.Robj2slave(y)
    mpi.bcast.cmd(output <- npdeptest(x, y, boot.num=2, method="summation"), caller.execute=TRUE)
    expect_s3_class(output, "deptest")
  })

  # 10. npsymtest
  test_that("npsymtest works in parallel", {
    n <- 50
    x <- rnorm(n)
    mpi.bcast.Robj2slave(x)
    mpi.bcast.cmd(output <- npsymtest(x, boot.num=2, method="summation"), caller.execute=TRUE)
    expect_s3_class(output, "symtest")
  })

  # 11. npunitest
  test_that("npunitest works in parallel", {
    n <- 50
    x <- rnorm(n)
    y <- rnorm(n)
    mpi.bcast.Robj2slave(x)
    mpi.bcast.Robj2slave(y)
    mpi.bcast.cmd(output <- npunitest(x, y, boot.num=2, method="summation"), caller.execute=TRUE)
    expect_s3_class(output, "unitest")
  })

  # 12. npcmstest / npqcmstest
  test_that("npcmstest works in parallel", {
    n <- 50
    x <- rnorm(n)
    y <- x + rnorm(n)
    mydat <- data.frame(x,y)
    mpi.bcast.Robj2slave(mydat)
    mpi.bcast.cmd(model <- lm(y~x, data=mydat, x=TRUE, y=TRUE), caller.execute=TRUE)
    mpi.bcast.cmd(output <- npcmstest(model = model, xdat = mydat$x, ydat = mydat$y, boot.num = 2), caller.execute=TRUE)
    expect_s3_class(output, "cmstest")
  })

  # 13. npconmode
  test_that("npconmode works in parallel", {
    data("Italy")
    Italy_sub <- Italy[1:50, ]
    mpi.bcast.Robj2slave(Italy_sub)
    mpi.bcast.cmd(bw <- npcdensbw(gdp~ordered(year), data=Italy_sub, bwmethod="normal-reference"), caller.execute=TRUE)
    mpi.bcast.cmd(model <- npconmode(bws=bw), caller.execute=TRUE)
    expect_s3_class(model, "conmode")
  })

  # 14. npuniden.reflect
  test_that("npuniden.reflect works in parallel", {
    n <- 50
    x <- sort(rbeta(n,5,1))
    mpi.bcast.Robj2slave(x)
    mpi.bcast.cmd(model <- npuniden.reflect(x), caller.execute=TRUE)
    expect_type(model, "list")
    expect_true("f" %in% names(model))
  })

})
