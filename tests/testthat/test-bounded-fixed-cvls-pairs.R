test_that("bounded fixed CVLS preserves both directions of resident pairs", {
  old <- options(np.messages=FALSE,np.tree=FALSE)
  on.exit(options(old),add=TRUE)
  pkg <- getNamespaceName(environment(npregbw));ev <- getFromNamespace(".npregbw_eval_only",pkg)
  set.seed(4040);n <- 60L;x <- data.frame(x=runif(n,-1,1));y <- sin(2*x$x)+rnorm(n,sd=.3);h <- .3
  for(cfg in list(c(0L,0L),c(1L,0L),c(2L,0L),c(2L,1L))) {
    degree <- cfg[1];bernstein <- as.logical(cfg[2]);H <- matrix(0,n,n)
    for(i in seq_len(n)) {
      z <- (x$x-x$x[i])/h;w <- dnorm(z);X <- outer(z,0:degree,"^")
      H[i,] <- solve(crossprod(X,w*X),t(X)*rep(w,each=ncol(X)))[1,]
    }
    fit <- drop(H%*%y);want <- mean(((y-fit)/(1-diag(H)))^2)
    for(bound in c("fixed","range")) for(tree in list(FALSE,TRUE,"auto")) {
      options(np.tree=tree)
      args <- list(xdat=x,ydat=y,bws=h,regtype=if(degree==0)"lc" else if(degree==1)"ll" else "lp",
        degree=degree,bernstein.basis=bernstein,ckerbound=bound,bwmethod="cv.ls",bandwidth.compute=FALSE)
      if(bound=="fixed"){args$ckerlb <- -1;args$ckerub <- 1}
      b <- do.call(npregbw,args)
      expect_lte(abs(ev(xdat=x,ydat=y,bws=b)$objective-want),2e-11)
      expect_equal(as.numeric(fitted(npreg(b,txdat=x,tydat=y,se=FALSE))),fit,tolerance=2e-11)
    }
  }
})

test_that("bounded pair correction composes with product and ordered kernels", {
  old <- options(np.messages=FALSE,np.tree=FALSE)
  on.exit(options(old),add=TRUE)
  pkg <- getNamespaceName(environment(npregbw));ev <- getFromNamespace(".npregbw_eval_only",pkg)
  set.seed(193);n <- 36L;xx <- runif(n,-1,1);y <- sin(2*xx)+rnorm(n,sd=.2)
  for(mixed in c(FALSE,TRUE)) {
    x <- if(mixed) data.frame(x=xx,z=ordered(rep(1:3,length.out=n))) else data.frame(x=xx,z=runif(n,-1,1))
    args <- list(xdat=x,ydat=y,bws=if(mixed)c(.35,.4) else c(.35,.5),regtype="lc",bwmethod="cv.ls",
      ckertype=if(mixed)"gaussian" else "epanechnikov",ckerbound="fixed",ckerlb=if(mixed)-1 else c(-1,-1),ckerub=if(mixed)1 else c(1,1),bandwidth.compute=FALSE)
    if(mixed)args$okertype <- "racineliyan"
    b <- do.call(npregbw,args)
    deleted <- vapply(seq_len(n),function(i)as.numeric(fitted(npreg(b,txdat=x[-i,,drop=FALSE],tydat=y[-i],exdat=x[i,,drop=FALSE],se=FALSE))),0.0)
    want <- mean((y-deleted)^2)
    for(tree in list(FALSE,TRUE,"auto")) {
      options(np.tree=tree)
      expect_lte(abs(ev(xdat=x,ydat=y,bws=b)$objective-want),2e-10)
    }
  }
})
