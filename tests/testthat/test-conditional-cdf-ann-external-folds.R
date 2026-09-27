test_that("external adaptive CDF objectives use literal deleted geometry", {
  old <- options(np.messages=FALSE,np.largeh=FALSE,np.largelambda=FALSE,
                 np.extendednn=FALSE)
  on.exit(options(old),add=TRUE)
  ns <- getNamespaceName(environment(npcdistbw))
  evaluate <- getFromNamespace(".npcdistbw_eval_only",ns)
  x <- data.frame(x=c(-.9,-.65,-.4,-.1,.15,.3,.55,.72,.94))
  y <- data.frame(y=c(-.5,-.3,.1,-.05,.3,.2,.65,.5,.82))
  grid <- data.frame(y=c(-.2,.2,.7,.2))
  n <- nrow(x); k <- 6L
  for(degree in 0:2) {
    loss <- 0
    for(i in seq_len(n)) {
      donors <- setdiff(seq_len(n),i)
      hx <- vapply(donors,function(j)
        sort(abs(x$x[-c(i,j)]-x$x[j]))[k],0.0)
      hy <- vapply(donors,function(j)
        sort(abs(y$y[-c(i,j)]-y$y[j]))[k],0.0)
      u <- (x$x[i]-x$x[donors])/hx
      w <- ifelse(abs(u)<sqrt(5),3/(4*sqrt(5))*(1-u*u/5)/hx,0)
      z <- vapply(0:degree,function(d)
        (x$x[donors]-x$x[i])^d,numeric(n-1L))
      influence <- w*drop(z%*%solve(crossprod(z,w*z),c(1,rep(0,degree))))
      for(q in grid$y) {
        fitted <- sum(influence*pnorm((q-y$y[donors])/hy))
        loss <- loss+((y$y[i]<=q)-fitted)^2
      }
    }
    expected <- loss/(n*nrow(grid))
    for(tree in c(FALSE,TRUE)) {
      options(np.tree=tree)
      args <- list(xdat=x,ydat=y,bws=c(k,k),bwtype="adaptive_nn",
        cxkertype="epanechnikov",cykertype="gaussian",
        regtype=if(degree)"lp"else"lc",bandwidth.compute=FALSE)
      if(degree) {args$degree<-degree;args$bernstein.basis<-TRUE}
      b <- do.call(npcdistbw,args)
      expect_equal(evaluate(x,y,b,gydat=grid,do.full.integral=FALSE)$objective,
                   expected,tolerance=2e-9)
    }
  }
})
