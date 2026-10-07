test_that("direct pooled single-index construction owns contrasts on the master", {
  skip_if(isTRUE(getOption("npRmpi.local.regression.mode", FALSE)),
          "This test requires distributed workers, not rank-local source mode")
  if (!spawn_mpi_slaves(1L)) skip("MPI pool unavailable")
  on.exit(close_mpi_slaves(force=TRUE),add=TRUE)
  old<-options(np.messages=FALSE,contrasts=c("contr.sum","contr.helmert"))
  on.exit(options(old),add=TRUE,after=FALSE)
  .r32.saved.contrasts <- mpi.remote.exec(getOption("contrasts"),simplify=FALSE)
  mpi.bcast.Robj2slave(.r32.saved.contrasts)
  on.exit(mpi.bcast.cmd({
    options(contrasts=.r32.saved.contrasts[[mpi.comm.rank(1L)]])
    rm(.r32.saved.contrasts)
  },caller.execute=FALSE),add=TRUE,after=FALSE)
  mpi.bcast.cmd(options(contrasts=c("contr.treatment","contr.poly")),caller.execute=FALSE)
  set.seed(619);n<-100L
  d<-data.frame(x=rnorm(n),x2=runif(n),g=factor(sample(c("A","B","C"),n,TRUE)))
  d$y<-sin(d$x+.5*d$x2+.4*(d$g=="B")-.6*(d$g=="C"))+rnorm(n,sd=.2)
  mm<-as.data.frame(model.matrix(y~x+x2+g,d)[,-1L,drop=FALSE])
  reference<-npindex(txdat=mm,tydat=d$y,nmulti=1L,se=FALSE)
  formula<-npindex(y~x+x2+g,data=d,nmulti=1L,se=FALSE)
  native<-npindex(txdat=d[c("x","x2","g")],tydat=d$y,nmulti=1L,se=FALSE)
  for(f in list(formula,native)) {
    expect_equal(unname(f$bws$beta),unname(reference$bws$beta),tolerance=1e-12)
    expect_equal(f$bw,reference$bw,tolerance=1e-12)
    expect_equal(fitted(f),fitted(reference),tolerance=1e-12)
    expect_equal(f$MSE,reference$MSE,tolerance=1e-12)
  }
  after<-mpi.remote.exec(getOption("contrasts"),simplify=FALSE)
  expect_true(all(vapply(after,function(x)identical(unname(x),c("contr.treatment","contr.poly")),logical(1))))
})
