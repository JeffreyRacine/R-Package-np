test_that("certified conditional bootstrap surfaces retain the point estimator", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force=TRUE),add=TRUE)
  old <- options(np.messages=FALSE)
  on.exit(options(old),add=TRUE)
  x <- data.frame(x=seq(-1,1,length.out=24))
  y <- data.frame(y=sin(seq_len(24)*.7))
  ex <- data.frame(x=rep(.13,7))
  ey <- data.frame(y=seq(-.6,.7,length.out=7))
  for(cdf in c(FALSE,TRUE)) for(type in c("lc","lp")) {
    bwfun <- if(cdf) npcdistbw else npcdensbw
    fitfun <- if(cdf) npcdist else npcdens
    b <- do.call(bwfun,c(list(xdat=x,ydat=y,bws=c(.4,.5),
      bandwidth.compute=FALSE,regtype=type),
      if(type=="lp") list(degree=0) else list()))
    for(gradient in c(FALSE,TRUE)) {
      args <- list(xdat=x,ydat=y,exdat=ex,eydat=ey,cdf=cdf,quantreg=FALSE,
        tau=.5,gradients=gradient,gradient.index=1L,slice.index=2L,
        plot.errors.boot.method="inid",plot.errors.boot.nonfixed="exact",
        plot.errors.boot.blocklen=NULL,plot.errors.boot.num=39L,
        plot.errors.center="estimate",plot.errors.type="pointwise",
        plot.errors.alpha=.05,bws=b)
      set.seed(274); raw <- do.call(compute.bootstrap.errors.conbandwidth,
                                  c(args,list(proper=FALSE)))
      seed <- .Random.seed
      set.seed(274); proper <- do.call(compute.bootstrap.errors.conbandwidth,
                                     c(args,list(proper=TRUE)))
      expect_identical(.Random.seed,seed)
      # MPI profiles measure two distinct executions, not estimator state.
      proper$timing.profile <- raw$timing.profile <- NULL
      expect_equal(proper,raw,tolerance=0)
      a <- fitfun(b,txdat=x,tydat=y,exdat=ex,eydat=ey,gradients=gradient,proper=FALSE)
      z <- fitfun(b,txdat=x,tydat=y,exdat=ex,eydat=ey,gradients=gradient,proper=TRUE)
      expect_equal(fitted(a),fitted(z),tolerance=0)
      if(gradient) expect_equal(gradients(a),gradients(z),tolerance=0)
    }
  }
})
