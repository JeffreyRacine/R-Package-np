test_that("conditional hybrid ratios reach the native objective owner", {
  old <- options(np.messages=FALSE,np.tree=FALSE,np.largeh=FALSE)
  on.exit(options(old),add=TRUE)
  x <- data.frame(x=c(.04,.11,.19,.3,.42,.56,.69,.78,.87,.96))
  y <- data.frame(y=c(.03,.08,.15,.23,.39,.52,.71,.82,.9,.98))
  build <- function(type="fixed",grid="hybrid",ratios=c(.2,.55,.25),yy=y) {
    npcdensbw(xdat=x,ydat=yy,
      bws=rep(if(type=="fixed") .16 else 4,ncol(yy)+1L),
      bwtype=type,bwmethod="cv.ls",bandwidth.compute=FALSE,
      cxkertype="beta",cxkerbound="fixed",cxkerlb=0,cxkerub=1,
      cykertype="beta",cykerbound="fixed",cykerlb=rep(0,ncol(yy)),cykerub=rep(1,ncol(yy)),
      cvls.quadrature.grid=grid,cvls.quadrature.points=c(12L,6L),
      cvls.quadrature.ratios=ratios)
  }
  value <- function(b,yy=y) .npcdensbw_eval_only(x,yy,b)$objective
  for(type in c("fixed","generalized_nn","adaptive_nn")) {
    expect_equal(value(build(type,ratios=c(1,0,0))),
                 value(build(type,grid="uniform")),tolerance=2e-12)
    for(grid in c("uniform","sample"))
      expect_identical(value(build(type,grid,c(1,0,0))),
                       value(build(type,grid,c(0,0,1))))
  }
  standard <- value(build())
  sample <- value(build(ratios=c(0,1,0)))
  expect_gt(abs(sample-standard),1e-6)
  expect_identical(value(build()),standard) # no prior-call global leakage
  yy <- data.frame(y=y$y,z=rev(y$y)^2)
  expect_error(build(yy=yy),"only for scalar continuous responses")
  expect_identical(value(build(grid="uniform",ratios=c(1,0,0),yy=yy),yy),
                   value(build(grid="uniform",ratios=c(0,0,1),yy=yy),yy))
  for(ratio in list(c(-.1,.6,.5),c(1,1,1),c(NA,0,1),c(Inf,0,0)))
    expect_error(build(ratios=ratio),"cvls.quadrature.ratios")
  # A short real search must consume the same uniform-only hybrid rule.
  set.seed(16)
  h <- npcdensbw(build(ratios=c(1,0,0)),xdat=x,ydat=y,
                 nmulti=1L,itmax=2L)
  set.seed(16)
  u <- npcdensbw(build(grid="uniform"),xdat=x,ydat=y,
                 nmulti=1L,itmax=2L)
  expect_equal(h$bw,u$bw,tolerance=2e-12)
  expect_equal(h$fval,u$fval,tolerance=2e-12)
})
