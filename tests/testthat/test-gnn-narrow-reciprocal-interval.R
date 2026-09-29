# Independent 70-digit Gaussian pair integrals on the NN radius partition;
# repeating at 100 digits gives the same reference. The ordinary fixture has
# distinct q endpoints whose double reciprocals coincide in the pair owner.
# The oracle integrates exp(-A*u^2/2-B*u-1)/(2*pi) in closed form on each piece.
test_that("GNN product overlaps retain narrow reciprocal intervals", {
  old <- options(np.messages=FALSE,np.tree=FALSE,np.extendednn=FALSE)
  on.exit(options(old),add=TRUE)
  skip_if_not(spawn_mpi_slaves(1L), "MPI session unavailable")
  on.exit(close_mpi_slaves(force=TRUE),add=TRUE)
  dat <- data.frame(y=c(-1.4,.2,.95,-.71,1.73,-.1,.57),
    y2=c(.42,-1.1,1.6,.03,-.63,.89,-.21),
    y3=c(-.38,1.42,-.04,-.91,.61,-1.73,.27))
  reference <- -0.092632597494878527119803973306912014907507880861878356854253414666
  for(tree in c(FALSE,TRUE)) {
    options(np.tree=tree)
    b <- npudensbw(dat=dat,bws=rep(2,3),bwtype="generalized_nn",
      bwmethod="cv.ls",bandwidth.compute=FALSE)
    command <- substitute({options(np.tree=TREE,np.messages=FALSE,np.extendednn=FALSE)
      get('npudensbw.bandwidth',asNamespace('npRmpi'))(dat=DAT,bws=B,
        eval.only=TRUE,invalid.penalty='dbmax',nmulti=1L)$fval},list(TREE=tree,DAT=dat,B=b))
    value <- get('.npRmpi_bcast_cmd_expr',asNamespace('npRmpi'))(
      command,comm=1L,caller.execute=TRUE)
    expect_true(is.finite(value) && abs(value-reference)<=1e-9)
  }
})
