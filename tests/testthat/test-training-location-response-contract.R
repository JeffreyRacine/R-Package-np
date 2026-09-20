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
                   tolerance=1e-11,info=paste(family,missing.row))
      args$eydat <- NULL
      unchanged <- do.call(get(family),args)
      expect_equal(fitted(actual),fitted(unchanged),tolerance=0)
    }
  }
})

test_that("explicit regression evaluation responses follow native tree row order", {
  old <- options(np.messages=FALSE,np.tree=getOption("np.tree"))
  on.exit(options(old),add=TRUE)
  set.seed(920317)
  n <- 35L
  x <- data.frame(x=runif(n),u=runif(n))
  y <- sin(5*x$x)+x$u+rnorm(n,sd=.1)
  for(type in c("fixed","generalized_nn","adaptive_nn")) {
    bw <- npregbw(xdat=x,ydat=y,bws=if(type=="fixed") c(.12,.16) else c(9,11),
                  bwtype=type,bandwidth.compute=FALSE)
    for(tree in list(FALSE,TRUE,"auto")) for(external in c(FALSE,TRUE)) {
      options(np.tree=tree)
      for(missing.row in c(FALSE,TRUE)) {
        xx <- x
        if(missing.row) xx$x[3] <- NA_real_
        keep <- complete.cases(xx)
        args <- list(bws=bw,txdat=xx,tydat=y,eydat=y+2,se=FALSE)
        if(external) args$exdat <- x
        fit <- do.call(npreg,args)
        response <- if(external) y+2 else y[keep]+2
        expect_equal(fit$MSE,mean((response-as.double(fitted(fit)))^2),
          tolerance=1e-11,info=paste(type,tree,external,missing.row))
        args$eydat <- NULL
        expect_equal(fitted(fit),fitted(do.call(npreg,args)),tolerance=0)
      }
    }
  }
})
