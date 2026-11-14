test_that("noncollective CDF degree recovery admits original-sample ordinary incumbents", {
  skip_if(.npRmpi_degree_search_collective(comm=1L), "requires the noncollective search owner")
  old<-options(np.messages=FALSE,np.tree=FALSE,np.extendednn=FALSE)
  on.exit(options(old),add=TRUE)
  recovered<-NULL
  local_mocked_bindings(.np_nomad_search=function(...) {
    recovered<<-list(...)$recover_start
    stop("R77 capture prepared recovery owner")
  },.package="npRmpi")
  x<-data.frame(x=c(rep(0,22),1,2));y<-data.frame(y=seq(.1,.9,length.out=24))
  for(full in c(FALSE,TRUE)) {
  expect_error(npcdistbw(xdat=x,ydat=y,do.full.integral=full,gydat=if(full)NULL else y[c(2,8,20),,drop=FALSE],
    bwtype="generalized_nn",regtype="lp",nomad=TRUE,degree.min=0L,
    degree.max=1L,degree.start=0L,degree.verify=FALSE,nomad.nmulti=2L,
    nmulti=1L,itmax=2L,nomad.opts=list(MAX_BB_EVAL=2L)),
    "R77 capture prepared recovery owner",fixed=TRUE)
  expect_true(is.function(recovered))
  out<-recovered(c(23,23,0))
  expect_true(out$found)
  expect_identical(out$point,c(22,22,0))
  expect_identical(out$evaluations,1L)
  expect_false(recovered(c(24,24,0))$found)
  raw<-mean(vapply(1:24,function(i) {
    d<-setdiff(1:24,i);hx<-sort(abs(x$x[d]-x$x[i]))[22]
    w<-dnorm((x$x[d]-x$x[i])/hx);w<-w/sum(w)
    mean(vapply(if(full)d else c(2,8,20),function(j) {
      hy<-sort(abs(y$y[d]-y$y[j]))[22]
      (as.numeric(y$y[i]<=y$y[j])-sum(w*pnorm((y$y[j]-y$y[d])/hy)))^2
    },0))
  },0))
  expect_lt(abs(out$objective-raw),2e-9)
  }
})
