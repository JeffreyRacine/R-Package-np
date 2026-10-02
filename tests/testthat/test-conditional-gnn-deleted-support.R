test_that('compact GNN CVLS preserves the adopted deleted-support policy', {
  old <- options(np.messages=FALSE,np.tree=FALSE)
  on.exit(options(old),add=TRUE)
  pkg <- getNamespaceName(environment(npcdensbw))
  evaluate <- function(x,y,b,tree) {
    command <- substitute(local({
      old <- options(np.tree=TREE,np.messages=FALSE)
      on.exit(options(old),add=TRUE)
      ev <- get('.npcdensbw_eval_only',asNamespace(PKG))
      args <- list(xdat=X,ydat=Y,bws=B,invalid.penalty='dbmax')
      if(PKG=='npRmpi')args$force.local <- FALSE
      do.call(ev,args)$objective
    }),list(X=data.frame(x=x),Y=y,B=b,TREE=tree,PKG=pkg))
    if(pkg=='npRmpi')get('.npRmpi_bcast_cmd_expr',asNamespace(pkg))(
      command,comm=1L,caller.execute=TRUE)else eval(command)
  }
  # n=100: uniform k=3 leaves only two nonzero donors after self deletion.
  # This is not restricted to a small training sample.
  set.seed(31);x <- runif(100,-1,1);y <- sin(2*x)+rnorm(100,sd=.4)
  # The incumbent compact-GNN admission requires at least degree+1
  # distinct nonzero-weight donors. A raw-coordinate ridge score here would
  # violate that policy and vary with predictor units (CF-212).
  for(bernstein in c(FALSE,TRUE)) {
    b <- npcdensbw(xdat=data.frame(x=x),ydat=y,bws=c(2,3),
      bwtype='generalized_nn',bwmethod='cv.ls',regtype='lp',degree=2L,
      bernstein.basis=bernstein,cxkertype='uniform',bandwidth.compute=FALSE)
    for(tree in list(FALSE,TRUE,'auto')) {
      expect_identical(evaluate(x,y,b,tree),-.Machine$double.xmax)
    }
  }
  # Retain the tied-data rejection control alongside the uniform support boundary.
  set.seed(30);x <- runif(30,-1,1);y <- sin(2*x)+rnorm(30,sd=.4)
  x[5] <- x[2];x[9] <- x[2];y[7] <- y[3];y[11] <- y[3]
  b <- npcdensbw(xdat=data.frame(x=x),ydat=y,bws=c(2,3),
    bwtype='generalized_nn',bwmethod='cv.ls',regtype='lp',degree=2L,
    cxkertype='epanechnikov',bandwidth.compute=FALSE)
  for(tree in list(FALSE,TRUE,'auto'))
    expect_identical(evaluate(x,y,b,tree),-.Machine$double.xmax)
  # Signed higher-order weights are donors, not invalid positive-mass tests.
  b <- npcdensbw(xdat=data.frame(x=x),ydat=y,bws=c(8,20),
    bwtype='generalized_nn',bwmethod='cv.ls',regtype='lp',degree=2L,
    cxkertype='epanechnikov',cxkerorder=8L,bandwidth.compute=FALSE)
  expect_lt(abs(evaluate(x,y,b,FALSE)),.Machine$double.xmax)
})
