.rly_matrix <- function(s, lambda) {
  a <- lambda^abs(outer(s,s,"-"))
  a/rowSums(a)
}

test_that("RLY operators use the actual retained support and donor normalizer", {
  old <- options(np.messages=FALSE,np.largeh=FALSE,np.largelambda=FALSE)
  on.exit(options(old))
  for(s in list(1:4,c(1,3,7),c(.5,2.5,6.5))) {
    ids <- rep(seq_along(s),seq_along(s)+2L)
    x <- data.frame(o=ordered(s[ids],levels=s))
    e <- data.frame(o=ordered(s,levels=s))
    for(lambda in c(0,.42,1-1e-9,1)) {
      K <- .rly_matrix(s,lambda)
      expected <- list(normal=K,convolution=tcrossprod(K),
                       integral=t(apply(K,1,cumsum)))
      for(op in names(expected)) {
        got <- npksum(txdat=x,exdat=e,bws=lambda,okertype="racineliyan",
                      operator=op,return.kernel.weights=TRUE)
        expect_equal(unname(got$kw),expected[[op]][ids,,drop=FALSE],
                     tolerance=2e-13,info=paste(op,lambda,paste(s,collapse=",")))
      }
      expect_equal(sum(fitted(npudens(tdat=x,edat=e,bws=lambda,
        okertype="racineliyan"))),1,tolerance=2e-13)
      expect_equal(as.numeric(fitted(npudist(tdat=x,edat=e,bws=lambda,
        okertype="racineliyan"))),cumsum(colMeans(K[ids,,drop=FALSE])),
        tolerance=2e-13)
    }
    lambda <- .42; eps <- 1e-6
    score <- npksum(txdat=x,exdat=e,bws=lambda,okertype="racineliyan",
                    compute.score=TRUE,return.kernel.weights=TRUE,
                    return.derivative.kernel.weights=TRUE)$p.kw
    reference <- (.rly_matrix(s,lambda+eps)-.rly_matrix(s,lambda-eps))/(2*eps)
    expect_equal(as.numeric(score),as.numeric(reference[ids,,drop=FALSE]),
                 tolerance=2e-9)
  }
})

test_that("RLY near-upper shortcuts respect the complete support", {
  old<-options(np.messages=FALSE,np.largeh=FALSE,np.largelambda=TRUE,
               np.disc.upper.rel.tol=.01)
  on.exit(options(old))
  s<-c(1,3,7);x<-data.frame(o=ordered(rep(s,each=5),levels=s))
  e<-data.frame(o=ordered(s,levels=s))
  for(lambda in c(.9899998,.99,.9900002,1)) {
    expected<-.rly_matrix(s,lambda)[rep(1:3,each=5),]
    got<-npksum(txdat=x,exdat=e,bws=lambda,okertype="racineliyan",
               return.kernel.weights=TRUE)$kw
    expect_equal(unname(got),expected,tolerance=2e-13)
  }
})

test_that("RLY density CV uses retained-support overlap and directed deletion", {
  old<-options(np.messages=FALSE,np.largeh=FALSE,np.largelambda=FALSE,
               np.categorical.compress=FALSE)
  on.exit(options(old))
  s<-c(1,3,7);ids<-rep(1:3,c(7,4,9));n<-length(ids)
  x<-data.frame(o=ordered(s[ids],levels=s));lambda<-.42
  K<-.rly_matrix(s,lambda)[ids,,drop=FALSE]
  fit<-colMeans(K); deleted<-K[,ids];diag(deleted)<-0
  loo<-colSums(deleted)/(n-1)
  for(compress in c(FALSE,TRUE))for(method in c("cv.ls","cv.ml")) {
    options(np.categorical.compress=compress)
    bw<-npudensbw(dat=x,bws=lambda,bwmethod=method,okertype="racineliyan",
                  bandwidth.compute=FALSE)
    got<-getFromNamespace("npudensbw.bandwidth","npRmpi")(dat=x,bws=bw,
      bandwidth.compute=TRUE,eval.only=TRUE,nmulti=1)$fval
    # fval stores the maximized score: -ISE for LS, total log LOO for ML.
    expected<-if(method=="cv.ls")2*mean(loo)-sum(fit^2) else sum(log(loo))
    expect_equal(as.numeric(got),expected,tolerance=2e-12)
  }
})

test_that("directed RLY paired CV and row-specific CVAIC match independent WLS", {
  old <- options(np.messages=FALSE,np.largeh=FALSE,np.largelambda=FALSE,np.tree=FALSE)
  on.exit(options(old))
  set.seed(620)
  n <- 140L; v <- sort(runif(n,-1,1)); s <- c(1,3,7)
  o <- sample(s,n,TRUE)
  x <- data.frame(v=v,o=ordered(o,levels=s))
  y <- sin(3*v)+o/5+rnorm(n,sd=.12)
  eval_cv <- getFromNamespace(".npregbw_eval_only","npRmpi")
  for(ker in c("gaussian","epanechnikov")) for(degree in 0:3) {
    h <- .35; lambda <- .42
    H <- t(vapply(seq_len(n),function(j) {
      z <- (v-v[j])/h
      w <- if(ker=="gaussian") dnorm(z) else pmax(1-z*z/5,0)*3/(4*sqrt(5))
      w <- w*lambda^abs(o-o[j])/vapply(o,function(a)sum(lambda^abs(a-s)),0)
      X <- outer(z,0:degree,"^")
      solve(crossprod(X,w*X),t(w*X))[1L,]
    },numeric(n)))
    fitted_ref <- drop(H%*%y); tr <- sum(diag(H))
    loo <- mean(((y-fitted_ref)/(1-diag(H)))^2)
    aic <- log(mean((y-fitted_ref)^2))+(1+tr/n)/(1-(tr+2)/n)
    for(tree in c(FALSE,TRUE)) for(method in c("cv.ls","cv.aic")) {
      options(np.tree=tree)
      bw <- npregbw(xdat=x,ydat=y,bws=c(h,lambda),bandwidth.compute=FALSE,
        regtype=if(degree==0)"lc" else "lp",degree=degree,degree.select="manual",
        ckertype=ker,okertype="racineliyan",bwmethod=method,bernstein.basis=FALSE)
      expect_equal(as.numeric(fitted(npreg(bws=bw,txdat=x,tydat=y))),fitted_ref,
                   tolerance=2e-10,info=paste(ker,degree,tree,method,"fit"))
      expect_equal(as.numeric(eval_cv(x,y,bw)$objective),
                   if(method=="cv.ls")loo else aic,tolerance=2e-9,
                   info=paste(ker,degree,tree,method,"objective"))
    }
  }
})
