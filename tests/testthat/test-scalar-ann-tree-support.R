test_that("scalar adaptive tree rows preserve literal deleted radii across tiles", {
  old <- options(np.messages=FALSE,np.largeh=FALSE,np.largelambda=FALSE,
                 np.extendednn=FALSE)
  on.exit(options(old),add=TRUE)
  ns <- getNamespaceName(environment(npregbw))
  evaluate <- getFromNamespace(".npregbw_eval_only",ns)
  set.seed(264)
  for(n in c(63L,64L,65L)) {
    x <- data.frame(x=runif(n,-1,1),z=runif(n,-1,1))
    y <- sin(x$x)+x$z
    k <- ceiling(n*.4)
    literal <- vapply(seq_len(n),function(i) {
      donors <- setdiff(seq_len(n),i)
      w <- rep(1,n-1L)
      for(d in seq_len(ncol(x))) {
        radius <- vapply(donors,function(j)
          sort(abs(x[-c(i,j),d]-x[j,d]))[k],0.0)
        u <- (x[i,d]-x[donors,d])/radius
        w <- w*ifelse(abs(u)<sqrt(5),
          3/(4*sqrt(5))*(1-u*u/5)/radius,0)
      }
      sum(w*y[donors])/sum(w)
    },0.0)
    expected <- mean((y-literal)^2)
    for(tree in list(FALSE,TRUE,"auto"))for(engine in c("lc","lp")) {
      options(np.tree=tree)
      args <- list(xdat=x,ydat=y,bws=rep(k,2),bwtype="adaptive_nn",
        regtype=engine,ckertype="epanechnikov",bandwidth.compute=FALSE)
      if(engine=="lp")args$degree <- c(0L,0L)
      b <- do.call(npregbw,args)
      expect_equal(evaluate(x,y,b)$objective,expected,tolerance=2e-12)
    }
  }
})
