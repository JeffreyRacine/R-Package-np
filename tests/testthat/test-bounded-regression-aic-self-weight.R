r23_aic_smoother <- function(x,kernel,order,type,bandwidth,degree,lower,upper) {
  n <- length(x)
  radius <- if(type=="fixed")rep(bandwidth,n) else vapply(seq_len(n),function(j)sort(abs(x[-j]-x[j]))[bandwidth],0.0)
  K <- function(z) switch(kernel,gaussian=dnorm(z)*if(order==4L)(1.5-.5*z^2) else 1,
    epanechnikov=ifelse(abs(z)<sqrt(5),3/(4*sqrt(5))*(1-z^2/5),0),uniform=.5*(abs(z)<1))
  F <- function(z) switch(kernel,gaussian=pnorm(z)+if(order==4L).5*z*dnorm(z) else 0,
    epanechnikov={v<-pmax(-sqrt(5),pmin(sqrt(5),z));.5+3/(4*sqrt(5))*(v-v^3/15)},uniform=pmax(0,pmin(1,(z+1)/2)))
  H <- matrix(0,n,n)
  for(i in seq_len(n)) {
    h <- if(type=="adaptive_nn")radius else radius[i]
    centre <- if(type=="adaptive_nn")x else x[i]
    mass <- F((upper-centre)/h)-F((lower-centre)/h)
    w <- K((x[i]-x)/h)/(h*mass)
    X <- outer(x-x[i],0:degree,"^")
    H[i,] <- solve(crossprod(X,w*X),t(X)*rep(w,each=ncol(X)))[1,]
  }
  H
}

test_that("bounded AIC self weights agree with their smoothing matrices", {
  skip_if_not(spawn_mpi_slaves(), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old <- options(np.messages=FALSE,np.tree=FALSE)
  on.exit(options(old),add=TRUE)
  pkg <- getNamespaceName(environment(npregbw));ev <- getFromNamespace(".npregbw_eval_only",pkg)
  set.seed(4040);n <- 60L;x <- data.frame(x=runif(n,-1,1));y <- sin(2*x$x)+rnorm(n,sd=.3)
  configs <- list(
    c("fixed","gaussian",0,2),c("fixed","gaussian",1,2),c("fixed","gaussian",2,2),
    c("fixed","epanechnikov",0,2),c("fixed","epanechnikov",2,2),c("fixed","uniform",0,2),
    c("generalized_nn","gaussian",0,2),c("generalized_nn","gaussian",2,2),c("generalized_nn","epanechnikov",0,2),
    c("adaptive_nn","gaussian",0,2),c("adaptive_nn","gaussian",2,2),c("adaptive_nn","epanechnikov",0,2),
    c("fixed","gaussian",0,4),c("adaptive_nn","gaussian",0,4))
  for(cfg in configs) {
    type <- cfg[1];kernel <- cfg[2];degree <- as.integer(cfg[3]);order <- as.integer(cfg[4])
    bw <- if(type=="fixed")if(order==4L).6 else .3 else if(order==4L)20L else 12L
    H <- r23_aic_smoother(x$x,kernel,order,type,bw,degree,-1,1)
    expected <- drop(H%*%y);tr <- sum(diag(H))
    want <- log(mean((y-expected)^2))+(1+tr/n)/(1-(tr+2)/n)
    for(tree in list(FALSE,TRUE,"auto")) {
      options(np.tree=tree)
      args <- list(xdat=x,ydat=y,bws=bw,bwtype=type,regtype=if(degree==0L)"lc" else if(degree==1L)"ll" else "lp",
        degree=degree,bernstein.basis=degree==2L,ckertype=kernel,ckerbound="fixed",ckerlb=-1,ckerub=1,
        bwmethod="cv.aic",bandwidth.compute=FALSE)
      if(kernel!="uniform")args$ckerorder <- order
      b <- do.call(npregbw,args)
      expect_equal(as.numeric(fitted(npreg(b,txdat=x,tydat=y,se=FALSE))),expected,tolerance=2e-11)
      expect_lte(abs(ev(xdat=x,ydat=y,bws=b)$objective-want),2e-11)
    }
  }
})

test_that("bounded AIC keeps retained range bounds and ordered self weights", {
  skip_if_not(spawn_mpi_slaves(), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old <- options(np.messages=FALSE,np.tree=FALSE)
  on.exit(options(old),add=TRUE)
  pkg <- getNamespaceName(environment(npregbw));ev <- getFromNamespace(".npregbw_eval_only",pkg)
  set.seed(191);n <- 48L;xx <- runif(n,-1,1);y <- sin(2*xx)+rnorm(n,sd=.2)
  for(mixed in c(FALSE,TRUE)) {
    x <- if(mixed)data.frame(x=xx,z=ordered(rep(1:3,length.out=n))) else data.frame(x=xx)
    args <- list(xdat=x,ydat=y,bws=if(mixed)c(.35,.4) else 12L,bwtype=if(mixed)"fixed" else "generalized_nn",
      regtype="lc",bwmethod="cv.aic",ckertype="gaussian",ckerbound="range",bandwidth.compute=FALSE)
    if(mixed)args$okertype <- "racineliyan"
    b <- do.call(npregbw,args)
    H <- if(mixed)npreghat(b,txdat=x,output="matrix") else r23_aic_smoother(xx,"gaussian",2L,"generalized_nn",12L,0L,min(xx),max(xx))
    expected <- drop(H%*%y);tr <- sum(diag(H))
    want <- log(mean((y-expected)^2))+(1+tr/n)/(1-(tr+2)/n)
    for(tree in list(FALSE,TRUE,"auto")) {
      options(np.tree=tree)
      expect_equal(as.numeric(fitted(npreg(b,txdat=x,tydat=y,se=FALSE))),expected,tolerance=2e-11)
      expect_lte(abs(ev(xdat=x,ydat=y,bws=b)$objective-want),2e-11)
    }
  }
})
