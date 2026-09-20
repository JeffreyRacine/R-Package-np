test_that("retained bandwidth routes reject recycling of paired rows", {
  old <- options(np.messages=FALSE)
  on.exit(options(old),add=TRUE)
  x <- data.frame(x=seq(.1,.9,length.out=20))
  y <- data.frame(y=sin(x$x))
  for (family in c("npcdens","npcdist")) {
    bw <- do.call(get(paste0(family,"bw")),list(xdat=x,ydat=y,bws=c(.2,.3),
      bandwidth.compute=FALSE))
    fun <- get(family)
    for (idx in list(1:10,rep(1:20,2))) {
      expect_error(fun(bw,txdat=x,tydat=y[idx,,drop=FALSE]),"same number of rows")
      expect_error(fun(bw,txdat=x,tydat=y,exdat=x,eydat=y[idx,,drop=FALSE]),
                   "same number of rows")
    }
    expect_length(fitted(fun(bw,txdat=x,tydat=y)),20L)
  }
  bw <- npscoefbw(xdat=x,ydat=y$y,zdat=x,bws=.3,bandwidth.compute=FALSE)
  for (idx in list(1:10,rep(1:20,2))) {
    expect_error(npscoef(bw,txdat=x,tydat=y$y[idx],tzdat=x,iterate=FALSE),
                 "same number of rows")
    expect_error(npscoef(bw,txdat=x,tydat=y$y,tzdat=x[idx,,drop=FALSE],iterate=FALSE),
                 "same number of rows")
    expect_error(npscoef(bw,txdat=x,tydat=y$y,tzdat=x,exdat=x,
      ezdat=x[idx,,drop=FALSE],iterate=FALSE),"same number of rows")
  }
  for (family in c("npindex","npplreg")) {
    args <- list(xdat=x,ydat=y$y,bandwidth.compute=FALSE,
      bws=if(family=="npindex") c(1,.3) else matrix(.3,2,1))
    if (family=="npplreg") args$zdat <- x
    bw <- do.call(get(paste0(family,"bw")),args)
    for (idx in list(1:10,rep(1:20,2))) {
      args <- list(bws=bw,txdat=x,tydat=y$y[idx],se=FALSE)
      if (family=="npplreg") args$tzdat <- x
      expect_error(do.call(get(family),args),"same number of rows")
    }
  }
})
