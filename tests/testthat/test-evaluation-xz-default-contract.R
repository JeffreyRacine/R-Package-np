test_that("separate evaluation X and Z roles honor their documented defaults", {
  old <- options(np.messages=FALSE)
  on.exit(options(old),add=TRUE)
  set.seed(920325)
  n <- 30L
  x <- data.frame(x=runif(n))
  z <- data.frame(z=runif(n))
  y <- x$x*(1+z$z)+sin(4*z$z)+rnorm(n,sd=.1)
  for (family in c("npscoef","npplreg")) {
    bw <- do.call(get(paste0(family,"bw")),list(xdat=x,zdat=z,ydat=y,
      bws=if(family=="npscoef") .3 else matrix(.3,2,1),bandwidth.compute=FALSE))
    for (omissions in c(FALSE,TRUE)) {
      xx <- x
      if (omissions) xx$x[3] <- NA_real_
      args <- list(bws=bw,txdat=xx,tydat=y,tzdat=z,se=FALSE)
      if(family=="npscoef") args$iterate <- FALSE
      fun <- get(family)
      ez <- data.frame(z=1-z$z)
      a <- do.call(fun,c(args,list(ezdat=ez)))
      b <- do.call(fun,c(args,list(exdat=xx,ezdat=ez)))
      expect_equal(fitted(a),fitted(b),tolerance=0)
      a <- do.call(fun,c(args,list(exdat=xx)))
      b <- do.call(fun,c(args,list(exdat=xx,ezdat=z)))
      expect_equal(fitted(a),fitted(b),tolerance=0)
      expect_error(do.call(fun,c(args,list(ezdat=ez[1:3,,drop=FALSE]))),
                   "same number of rows")
    }
  }
  bw <- npscoefbw(xdat=x,ydat=y,bws=.3,bandwidth.compute=FALSE)
  ez <- data.frame(x=1-x$x)
  a <- npscoef(bw,txdat=x,tydat=y,ezdat=ez,iterate=FALSE)
  b <- npscoef(bw,txdat=x,tydat=y,exdat=x,ezdat=ez,iterate=FALSE)
  expect_equal(fitted(a),fitted(b),tolerance=0)
  # Explicit-Z construction is an independent representation of the same model.
  bz <- npscoefbw(xdat=x,zdat=x,ydat=y,bws=.3,bandwidth.compute=FALSE)
  c <- npscoef(bz,txdat=x,tzdat=x,tydat=y,exdat=x,ezdat=ez,iterate=FALSE)
  expect_equal(fitted(a),fitted(c),tolerance=1e-11)
})
