test_that("hybrid conditional quadrature keeps support endpoints inside", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force=TRUE),add=TRUE)
  old <- options(np.messages=FALSE,np.tree=FALSE,np.largeh=FALSE)
  on.exit(options(old),add=TRUE)
  set.seed(920313)
  x <- data.frame(x=5*runif(18,.03,.97)+7)
  y <- data.frame(y=.9*runif(18,.03,.97)+.1)
  objective <- function(ub,q,kernel,grid="hybrid",ratios=c(.2,.55,.25)) {
    bw <- npcdensbw(xdat=x,ydat=y,bws=c(.18,1.5),
      bandwidth.compute=FALSE,bwmethod="cv.ls",regtype="lc",
      cxkertype="beta",cxkerbound="fixed",cxkerlb=7,cxkerub=12,
      cykertype=kernel,cykerbound="fixed",cykerlb=.1,cykerub=ub,
      cvls.quadrature.grid=grid,cvls.quadrature.points=c(q,15L),
      cvls.quadrature.ratios=ratios)
    .npcdensbw_eval_only(x,y,bw)$objective
  }
  for(kernel in c("beta","gaussian")) for(q in c(8L,15L,40L,75L,145L)) {
    a <- objective(1,q,kernel)
    b <- objective(1+4*.Machine$double.eps,q,kernel)
    expect_equal(a,b,tolerance=2e-12)
  }
})
