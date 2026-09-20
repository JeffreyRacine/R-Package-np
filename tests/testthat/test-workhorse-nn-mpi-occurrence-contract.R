test_that("MPI numerical NN rows preserve the corrected local objective", {
  skip_if_not(spawn_mpi_slaves(), "MPI pool unavailable")
  on.exit(close_mpi_slaves(), add = TRUE)
  withr::local_options(list(np.messages=FALSE,np.largeh=FALSE,np.largelambda=FALSE))
  set.seed(920499)
  x <- data.frame(x=runif(32,.02,.98))
  y <- sin(6*x$x)+rnorm(32,sd=.2)
  yy <- data.frame(y=(y-min(y)+.1)/(max(y)-min(y)+.2))
  for(type in c("generalized_nn","adaptive_nn")) {
    for(reg in c("lc","ll","lp")) for(criterion in c("cv.ls","cv.aic")) {
      args <- list(xdat=x,ydat=y,bws=10,bwtype=type,ckertype="beta",
        ckerbound="fixed",ckerlb=0,ckerub=1,regtype=reg,bwmethod=criterion,
        bandwidth.compute=FALSE)
      if(reg=="lp") args$degree <- 2L
      b <- do.call(npregbw,args)
      local <- .npregbw_eval_only(x,y,b)$objective
      cmd <- substitute(npRmpi:::.npregbw_eval_only(X,Y,B,localize=FALSE),
                        list(X=x,Y=y,B=b))
      distributed <- do.call(mpi.bcast.cmd,
        list(cmd=cmd,comm=1L,caller.execute=TRUE))$objective
      expect_equal(distributed,local,tolerance=1e-11,
                   info=paste(type,reg,criterion))
    }
    for(reg in c("lc","ll","lp")) for(kernels in list(
      c("beta","beta"),c("beta","gaussian"),c("gaussian","beta"))) {
      args <- list(xdat=x,ydat=yy,bws=c(10,10),bwtype=type,regtype=reg,
        cxkertype=kernels[1],cykertype=kernels[2],bwmethod="cv.ml",
        bandwidth.compute=FALSE)
      if(reg=="lp") args$degree <- 2L
      if(kernels[1]=="beta") args <- c(args,list(cxkerbound="fixed",cxkerlb=0,cxkerub=1))
      if(kernels[2]=="beta") args <- c(args,list(cykerbound="fixed",cykerlb=0,cykerub=1))
      b <- do.call(npcdensbw,args)
      local <- .npcdensbw_eval_only(x,yy,b)$objective
      cmd <- substitute(npRmpi:::.npcdensbw_eval_only(X,Y,B,force.local=FALSE),
                        list(X=x,Y=yy,B=b))
      distributed <- do.call(mpi.bcast.cmd,
        list(cmd=cmd,comm=1L,caller.execute=TRUE))$objective
      expect_equal(distributed,local,tolerance=1e-11,info=paste(type,reg,kernels))
    }
  }
})
