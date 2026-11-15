test_that('U5 receives the public plotting seed for fitted and bandwidth routes', {
  if (!spawn_mpi_slaves()) skip('MPI pool unavailable')
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  pkg <- getNamespaceName(environment(.np_udens_u5_prepare))
  withr::local_options(np.messages=FALSE, np.plot.progress=FALSE, np.tree=FALSE)
  original <- .np_udens_u5_prepare
  seen <- integer()
  testthat::local_mocked_bindings(.np_udens_u5_prepare=function(bws,x,random.seed=42L) {
    seen <<- c(seen,random.seed)
    original(bws,x,random.seed)
  }, .package=pkg)
  x <- data.frame(x=qbeta((seq_len(40L)-.5)/40L,2,3))
  for(order in c(4L,6L,8L)) {
    bw <- npudensbw(dat=x,bws=.12,bandwidth.compute=FALSE,
      ckertype='epanechnikov',ckerorder=order,ckerbound='fixed',ckerlb=0,ckerub=1)
    fit <- npudens(bws=bw,tdat=x,se=FALSE)
    for(object in list(fit,bw)) {
      seen <- integer()
      saved <- serialize(object,NULL)
      set.seed(82); before <- .Random.seed
      out <- plot(object,xdat=x,errors='asymptotic',center='bias-corrected',
        band='pointwise',random.seed=43L,neval=3L,output='data')
      expect_identical(seen,43L)
      expect_identical(out$request$random_seed,43L)
      expect_identical(before,.Random.seed)
      expect_identical(saved,serialize(object,NULL))
    }
  }
})
