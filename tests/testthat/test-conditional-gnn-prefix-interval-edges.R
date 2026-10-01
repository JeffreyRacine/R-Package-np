test_that('conditional GNN prefix retains finite response interval edges', {
  old <- options(np.messages=FALSE, np.macMseries.accelerate=FALSE, np.tree=FALSE)
  on.exit(options(old), add=TRUE)
  x <- data.frame(x=c(-.83,.66,-.21,.32,-.49,.93,.09))
  # Minimal regression fixtures, not claims about practical incidence. Targets
  # are independent deleted-sample WLS/full-line integrals, not np outputs.
  cases <- list(
    list(y=c(-.38,1.42,-.04,-.91,.61,-1.73,.27), k=2, kernel='gaussian',
         target=.059801807063833645),
    list(y=c(-1.4,.2,.95,-.71,1.73,-.1,.57), k=1, kernel='uniform', target=0),
    list(y=c(-1.4,.2,.2,.2,1.73,-.1,.57), k=3, kernel='uniform',
         target=.35435239719726486))
  ev <- get('.npcdensbw_eval_only', asNamespace('np'))
  for(z in cases) for(tree in c(FALSE,TRUE)) {
    options(np.tree=tree)
    b <- npcdensbw(xdat=x,ydat=z$y,bws=c(z$k,4),bwtype='generalized_nn',
      bwmethod='cv.ls',regtype='lc',cxkertype='epanechnikov',
      cykertype=z$kernel,bandwidth.compute=FALSE)
    value <- ev(x,z$y,b,invalid.penalty='dbmax')$objective
    expect_true(is.finite(value))
    expect_lte(abs(value-z$target),1e-9)
  }
  # An infinite reciprocal limit is not in general a finite integral.
  b <- npcdensbw(xdat=x,ydat=cases[[2]]$y,bws=c(1,4),bwtype='generalized_nn',
    bwmethod='cv.ls',regtype='lc',cxkertype='epanechnikov',
    cykertype='gaussian',bandwidth.compute=FALSE)
  expect_identical(ev(x,cases[[2]]$y,b,invalid.penalty='dbmax')$objective,
                   -.Machine$double.xmax)
})
