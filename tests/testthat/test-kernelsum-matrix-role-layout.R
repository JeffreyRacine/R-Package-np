test_that("kernel sums use weight-by-response layout across ordinary owners", {
  old <- options(np.messages=FALSE, np.largeh=FALSE, np.largelambda=FALSE)
  on.exit(options(old), add=TRUE)
  set.seed(30921)
  x <- data.frame(x=runif(23)); e <- data.frame(x=c(.11,.33,.66,.89))
  Y <- cbind(sin(x$x), 1+x$x, x$x^2-.3)
  W <- cbind(cos(3*x$x), x$x+.2)
  for (tree in c(FALSE, TRUE)) for (type in c("fixed","generalized_nn","adaptive_nn"))
    for (kernel in c("gaussian","epanechnikov","uniform")) {
      options(np.tree=tree)
      a <- list(txdat=x, exdat=e, bws=if(type=="fixed") .27 else 8,
                bwtype=type, ckertype=kernel)
      K <- do.call(npksum, c(a,list(return.kernel.weights=TRUE)))$kw
      for (power in 1:2) {
        got <- do.call(npksum,c(a,list(tydat=Y,weights=W,kernel.pow=power)))$ksum
        want <- array(0,c(ncol(W),ncol(Y),nrow(e)))
        for (j in seq_len(nrow(e))) want[,,j] <- crossprod(W,Y*K[,j]^power)
        expect_equal(got,want,tolerance=2e-12)
      }
    }
  beta <- npksum(txdat=x,exdat=e,tydat=Y,weights=W,bws=.08,ckertype="beta",
                ckerbound="fixed",ckerlb=0,ckerub=1,return.kernel.weights=TRUE)
  want <- array(0,c(2L,3L,nrow(e)))
  for (j in seq_len(nrow(e))) want[,,j] <- crossprod(W,Y*beta$kw[,j])
  expect_equal(beta$ksum,want,tolerance=2e-12)
})

test_that("matrix role decoding covers permutations and paired powers", {
  old <- options(np.messages=FALSE,np.tree=TRUE,np.largeh=FALSE,np.largelambda=FALSE)
  on.exit(options(old),add=TRUE)
  set.seed(34921); n<-24L
  x<-data.frame(z=runif(n),o=ordered(rep(0:2,length.out=n)),
                u=factor(rep(c("a","b"),length.out=n)))
  Y<-cbind(1+x$z, sin(3*x$z), x$z^2); W<-cbind(cos(x$z),seq_len(n)/n)
  e<-x[c(2,7,12),,drop=FALSE]; h<-c(.32,.4,.2)
  for (kind in c("derivative","score","ocg")) {
    a<-list(txdat=x,exdat=e,bws=h,tydat=Y,weights=W)
    if(kind=="derivative")a$permutation.operator<-"derivative"
    if(kind=="score")a$compute.score<-TRUE
    if(kind=="ocg")a$compute.ocg<-TRUE
    got<-do.call(npksum,a)
    q<-if(kind=="derivative")1L else 2L
    want<-array(0,c(2L,3L,3L,q))
    for(i in 1:2)for(j in 1:3) {
      scalar<-a;scalar$tydat<-Y[,j];scalar$weights<-W[,i,drop=FALSE]
      want[i,j,,]<-do.call(npksum,scalar)$p.ksum
    }
    dim(want)<-dim(want)[dim(want)>1L]
    expect_equal(got$p.ksum,want,tolerance=2e-12)
  }
  a<-list(txdat=x,exdat=e,bws=h)
  both<-.npksum_power12(txdat=x,exdat=e,bws=h)
  expect_identical(both$ksum,do.call(npksum,a)$ksum)
  expect_identical(both$ksum.power2,do.call(npksum,c(a,list(kernel.pow=2)))$ksum)
  for(type in c("fixed","generalized_nn","adaptive_nn")) {
    bw<-c(if(type=="fixed").32 else 7,.4,.2)
    bw<-kbandwidth.numeric(bw=bw,bwtype=type,nobs=n,xdati=untangle(x),xnames=names(x))
    got<-.np_estimator_loo_ksum(txdat=x,tydat=Y,weights=W,bws=bw,
                               bwtype=type,leave.one.out=TRUE)$ksum
    want<-array(0,c(2L,3L,n))
    for(i in 1:2)for(j in 1:3)
      want[i,j,]<-.np_estimator_loo_ksum(txdat=x,tydat=Y[,j],weights=W[,i,drop=FALSE],
                            bws=bw,bwtype=type,leave.one.out=TRUE)$ksum
    expect_equal(got,want,tolerance=2e-12)
  }
})

test_that("smooth coefficient backfitting keeps its scalar update contract", {
  old<-options(np.messages=FALSE,np.tree=FALSE,np.largeh=FALSE,np.largelambda=FALSE)
  on.exit(options(old),add=TRUE)
  set.seed(882);n<-24L;x<-data.frame(x=runif(n,.2,1.4));z<-data.frame(z=runif(n,.1,.9))
  y<-1+z$z+sin(z$z*3)*x$x+rnorm(n,sd=.06)
  bw<-suppressWarnings(npscoefbw(xdat=x,zdat=z,ydat=y,bws=.25,
    bandwidth.compute=TRUE,cv.iterate=TRUE,backfit.iterate=TRUE,
    cv.num.iterations=1L,backfit.maxiter=2L,optim.maxit=2L,nmulti=1L,random.seed=11L))
  expect_false(is.null(bw$bw.fitted))
  initial<-npscoef(bws=bw,txdat=x,tzdat=z,tydat=y,iterate=FALSE,betas=TRUE)
  actual<-suppressWarnings(npscoef(bws=bw,txdat=x,tzdat=z,tydat=y,iterate=TRUE,
                                  maxiter=1L,tol=1e-12,betas=TRUE))
  beta<-coef(initial);W<-cbind(1,x$x);resid<-y-rowSums(W*beta)
  for(j in seq_len(ncol(W))) {
    partial<-W[,j]*beta[,j]+resid
    K<-npksum(txdat=z,bws=as.numeric(bw$bw.fitted[,j]),return.kernel.weights=TRUE)$kw
    beta[,j]<-as.numeric(crossprod(K,partial*W[,j]))/colSums(K*W[,j]^2)
    resid<-partial-W[,j]*beta[,j]
  }
  expect_equal(unname(coef(actual)),unname(beta),tolerance=1e-11)
  expect_equal(as.numeric(fitted(actual)),rowSums(W*beta),tolerance=1e-11)
})

test_that("IV local moment consumers follow W-by-Y kernel sums", {
  old<-options(np.messages=FALSE,np.tree=FALSE,np.largeh=FALSE,np.largelambda=FALSE)
  on.exit(options(old),add=TRUE)
  # Exercise the actual nested owners without running a full IV bandwidth search.
  owner<-get("npregiv.default",envir=environment(npregiv))
  env<-new.env(parent=environment(owner));env$p<-2L
  for(expr in as.list(body(owner))[-1L])
    if(is.call(expr)&&identical(expr[[1]],as.name("<-"))&&is.symbol(expr[[2]])&&
       as.character(expr[[2]]) %in% c("glpreg","minimand.cv.ls","minimand.cv.aic"))
      eval(expr,env)
  x<-data.frame(x=seq(-.9,.9,length.out=28));y<-sin(x$x*2)+.2*x$x^3
  W<-W.lp(xdat=x,degree=2L);h<-.65
  oracle<-function(loo) {
    K<-outer(x$x,x$x,function(a,b)dnorm((a-b)/h)/h)
    if(loo)diag(K)<-0
    H<-matrix(0,nrow(x),nrow(x))
    for(j in seq_len(nrow(x)))
      H[j,]<-as.numeric(W[j,,drop=FALSE]%*%solve(crossprod(W,W*K[,j]),t(W)*rep(K[,j],each=ncol(W))))
    list(mean=as.numeric(H%*%y),H=H)
  }
  for(loo in c(FALSE,TRUE)) {
    want<-oracle(loo)
    got<-env$glpreg(tydat=y,txdat=x,bws=h,degree=2L,leave.one.out=loo)
    expect_equal(got$mean,want$mean,tolerance=1e-10)
  }
  want<-oracle(TRUE)
  expect_equal(env$minimand.cv.ls(bws=h,ydat=y,xdat=x,degree=2L,W=W),
               mean((y-want$mean)^2),tolerance=1e-10)
  want<-oracle(FALSE);trH<-sum(diag(want$H));n<-nrow(x)
  expect_equal(env$minimand.cv.aic(bws=h,ydat=y,xdat=x,degree=2L,W=W),
    log(mean((y-want$mean)^2))+(1+trH/n)/(1-(trH+2)/n),tolerance=1e-10)
})
