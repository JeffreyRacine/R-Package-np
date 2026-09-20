test_that("explicit responses at training locations retain row ownership", {
  old <- options(np.messages=FALSE)
  on.exit(options(old),add=TRUE)
  set.seed(920317)
  n <- 35L
  x <- data.frame(x=runif(n),u=runif(n))
  z <- data.frame(z=runif(n))
  y <- sin(5*z$z)+x$x+rnorm(n,sd=.1)
  for (family in c("npreg","npindex","npscoef","npplreg")) {
    args <- list(xdat=x,ydat=y,bandwidth.compute=FALSE,
      bws=switch(family,npreg=c(.3,.4),npindex=c(1,.5,.3),
        npscoef=.3,npplreg=matrix(.3,3,1)))
    if (family %in% c("npscoef","npplreg")) args$zdat <- z
    bw <- do.call(get(paste0(family,"bw")),args)
    for (missing.row in c(FALSE,TRUE)) {
      xx <- x
      if (missing.row) xx$x[3] <- NA_real_
      keep <- complete.cases(xx,y,z)
      args <- list(bws=bw,txdat=xx,tydat=y,eydat=y+2,se=FALSE)
      if (family %in% c("npscoef","npplreg")) args$tzdat <- z
      if (family=="npscoef") args$iterate <- FALSE
      actual <- do.call(get(family),args)
      expect_equal(actual$MSE,mean((y[keep]+2-as.double(fitted(actual)))^2),
                   tolerance=1e-11)
      args$eydat <- NULL
      unchanged <- do.call(get(family),args)
      expect_equal(fitted(actual),fitted(unchanged),tolerance=0)
    }
  }
})
