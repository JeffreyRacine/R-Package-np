r23_ann_weights <- function(train, evaluation, k, kernel, lower, upper) {
  h <- vapply(seq_along(train),function(j) sort(abs(train[-j]-train[j]))[k],0.0)
  u <- outer(evaluation,train,"-")/rep(h,each=length(evaluation))
  cdf <- switch(kernel,
    gaussian=pnorm,
    epanechnikov=function(z) {v<-pmax(-sqrt(5),pmin(sqrt(5),z)); .5+3/(4*sqrt(5))*(v-v^3/15)},
    uniform=function(z) pmax(0,pmin(1,(z+1)/2)))
  K <- switch(kernel,gaussian=dnorm(u),
    epanechnikov=ifelse(abs(u)<sqrt(5),3/(4*sqrt(5))*(1-u^2/5),0),
    uniform=.5*(abs(u)<1))
  mass <- cdf((upper-train)/h)-cdf((lower-train)/h)
  sweep(K,2L,h*mass,"/")
}

test_that("bounded ANN scalar objectives retain the fitted donor convention", {
  skip_if_not(spawn_mpi_slaves(), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old <- options(np.messages=FALSE,np.tree=FALSE)
  on.exit(options(old),add=TRUE)
  pkg <- getNamespaceName(environment(npregbw))
  evaluate <- getFromNamespace(".npregbw_eval_only",pkg)
  set.seed(4040); n <- 36L
  x <- data.frame(x=runif(n,-1,1));y <- sin(2*x$x)+rnorm(n,sd=.3);k <- 9L
  for(kernel in c("gaussian","epanechnikov","uniform")) for(bound in c("fixed","range")) {
    options(np.tree=FALSE)
    args <- list(xdat=x,ydat=y,bws=k,bwtype="adaptive_nn",regtype="lc",bwmethod="cv.ls",
      ckertype=kernel,ckerbound=bound,bandwidth.compute=FALSE)
    lower <- if(bound=="fixed") -1 else min(x$x)
    upper <- if(bound=="fixed") 1 else max(x$x)
    if(bound=="fixed") {args$ckerlb<-lower;args$ckerub<-upper}
    b <- do.call(npregbw,args)
    deleted <- vapply(seq_len(n),function(i) {
      w <- r23_ann_weights(x$x[-i],x$x[i],k,kernel,lower,upper)
      sum(w*y[-i])/sum(w)
    },0.0)
    for(i in c(1L,18L,36L)) {
      f <- fitted(npreg(b,txdat=x[-i,,drop=FALSE],tydat=y[-i],exdat=x[i,,drop=FALSE],se=FALSE))
      expect_equal(as.numeric(f),deleted[i],tolerance=2e-11)
    }
    W <- r23_ann_weights(x$x,x$x,k,kernel,lower,upper)
    expect_equal(as.numeric(fitted(npreg(b,txdat=x,tydat=y,se=FALSE))),drop(W%*%y)/rowSums(W),tolerance=2e-11)
    want <- mean((y-deleted)^2)
    for(tree in list(FALSE,TRUE,"auto")) {
      options(np.tree=tree)
      got <- evaluate(xdat=x,ydat=y,bws=b)$objective
      expect_true(is.finite(got))
      expect_lte(abs(got-want),2e-11)
    }
  }
})
