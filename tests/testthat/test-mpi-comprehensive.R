context("MPI Comprehensive Examples")

test_that("Comprehensive MPI examples work correctly", {
  # # skip_on_cran()
  
  if (!spawn_mpi_slaves()) {
    skip("Could not spawn MPI slaves for testing")
  }

  # Ensure cleanup on exit
  on.exit(try(close_mpi_slaves(force=TRUE), silent=TRUE))
  options(npRmpi.autodispatch = TRUE, np.messages = FALSE)

  # 1. npudensbw / npudens (Unconditional density)
  test_that("npudens works in parallel", {
    n <- 50
    x <- rnorm(n)
    bw <- npudensbw(~x, bwmethod="normal-reference")
    fhat <- npudens(bw)
    expect_s3_class(fhat, "npdensity")
  })

  # 2. npcdensbw / npcdens (Conditional density)
  test_that("npcdens works in parallel", {
    data("Italy")
    Italy_sub <- Italy[1:50, ]
    bw <- npcdensbw(gdp~ordered(year), data=Italy_sub, bwmethod="normal-reference")
    fhat <- npcdens(bw)
    expect_s3_class(fhat, "condensity")
  })

  # 3. npudistbw / npudist (Unconditional distribution)
  test_that("npudist works in parallel", {
    n <- 50
    x <- rnorm(n)
    bw <- npudistbw(~x, bwmethod="normal-reference")
    Fhat <- npudist(bw)
    expect_s3_class(Fhat, "npdistribution")
  })

  # 4. npcdistbw / npcdist (Conditional distribution)
  test_that("npcdist works in parallel", {
    data("Italy")
    Italy_sub <- Italy[1:50, ]
    bw <- npcdistbw(gdp~ordered(year), data=Italy_sub, bwmethod="normal-reference")
    Fhat <- npcdist(bw)
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
    # bwmethod normal-reference not allowed for npplregbw (it uses npregbw internally which might but scbandwidth doesn't)
    # npregbw default is cv.ls
    bw <- npplregbw(formula=y~x1+x2|z1+z2, data=mydat, bandwidth.compute=FALSE, bws=matrix(0.1, 5, 2))
    pl <- npplreg(bws=bw)
    expect_s3_class(pl, "plregression")
  })

  # 6. npindexbw / npindex (Single index model)
  test_that("npindex works in parallel", {
    n <- 50
    x1 <- runif(n, min=-1, max=1)
    x2 <- runif(n, min=-1, max=1)
    y <- x1 - x2 + rnorm(n)
    mydat <- data.frame(x1,x2,y)
    # Use small nmulti for speed
    bw <- npindexbw(formula=y~x1+x2, data=mydat, nmulti=1, bandwidth.compute=FALSE, bws=c(1, 1, 0.1))
    model <- npindex(bws=bw)
    expect_s3_class(model, "singleindex")
  })

  # 7. npscoefbw / npscoef (Smooth coefficient)
  test_that("npscoef works in parallel", {
    n <- 50
    x <- runif(n)
    z <- runif(n, min=-2, max=2)
    y <- x*exp(z)*(1.0+rnorm(n,sd = 0.2))
    mydat <- data.frame(x,y,z)
    # bwmethod="normal-reference" not allowed for scbandwidth
    bw <- npscoefbw(y~x|z, data=mydat, bandwidth.compute=FALSE, bws=0.1)
    model <- npscoef(bws=bw)
    expect_s3_class(model, "smoothcoefficient")
  })

  # 8. npsigtest (Significance test)
  test_that("npsigtest works in parallel", {
    n <- 50
    z <- factor(rbinom(n,1,.5))
    x1 <- rnorm(n)
    y <- x1 + rnorm(n)
    mydat <- data.frame(z,x1,y)
    model <- npreg(y~z+x1, data=mydat, bws=c(0.1, 0.1), bandwidth.compute=FALSE)
    output <- npsigtest(model, boot.num=9)
    expect_s3_class(output, "sigtest")
  })

  # 9. npdeptest / npsdeptest
  test_that("npdeptest works in parallel", {
    n <- 50
    x <- rnorm(n)
    y <- x + rnorm(n)
    output <- npdeptest(x, y, boot.num=9, method="summation")
    expect_s3_class(output, "deptest")
  })

  # 10. npsymtest
  test_that("npsymtest works in parallel", {
    n <- 50
    x <- rnorm(n)
    output <- npsymtest(x, boot.num=9, method="summation")
    expect_s3_class(output, "symtest")
  })

  # 11. npunitest
  test_that("npunitest works in parallel", {
    n <- 50
    x <- rnorm(n)
    y <- rnorm(n)
    output <- npunitest(x, y, boot.num=9, method="summation")
    expect_s3_class(output, "unitest")
  })

  # 12. npcmstest / npqcmstest
  test_that("npcmstest works in parallel", {
    n <- 50
    x <- rnorm(n)
    y <- x + rnorm(n)
    mydat <- data.frame(x,y)
    model <- lm(y~x, data=mydat, x=TRUE, y=TRUE)
    output <- npcmstest(model = model, xdat = mydat$x, ydat = mydat$y, boot.num = 9)
    expect_s3_class(output, "cmstest")
  })

  # 13. npconmode
  test_that("npconmode works in parallel", {
    data("Italy")
    Italy_sub <- Italy[1:50, ]
    # npconmode requires TYDAT to be discrete
    Italy_sub$gdp_cat <- factor(ifelse(Italy_sub$gdp > median(Italy_sub$gdp), "high", "low"))
    bw <- npcdensbw(gdp_cat~ordered(year), data=Italy_sub, bwmethod="normal-reference")
    model <- npconmode(bws=bw)
    expect_s3_class(model, "conmode")
  })

  # 14. npuniden.reflect
  test_that("npuniden.reflect works in parallel", {
    n <- 50
    x <- sort(rbeta(n,5,1))
    model <- npuniden.reflect(x)
    expect_type(model, "list")
    expect_true("f" %in% names(model))
  })

})
