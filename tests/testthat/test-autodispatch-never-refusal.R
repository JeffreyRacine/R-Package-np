test_that("option agreement refusal preserves existing mode normalization", {
  ns<-asNamespace('npRmpi');env<-environment(get('.npRmpi_distributed_call_impl',ns))
  check<-get('.npRmpi_autodispatch_require_option_sync',env);mode<-get('.npRmpi_autodispatch_option_mode',env)
  old<-options();on.exit(options(old),add=TRUE)
  for(value in list('onchange','always','on','al',NULL,NA_character_,1,c('never','always'))) {
    options(npRmpi.autodispatch.option.sync=value);expect_null(check())
  }
  for(value in c('never','ne')) {
    options(npRmpi.autodispatch.option.sync=value);seed<-if(exists('.Random.seed',.GlobalEnv))get('.Random.seed',.GlobalEnv)else NULL
    snapshot<-getOption('npRmpi.autodispatch.option.snapshot');expect_identical(mode(),'never')
    expect_error(check(),'automatic numerical calls require verified option agreement',fixed=TRUE)
    expect_identical(getOption('npRmpi.autodispatch.option.sync'),value)
    expect_identical(getOption('npRmpi.autodispatch.option.snapshot'),snapshot)
    expect_identical(if(exists('.Random.seed',.GlobalEnv))get('.Random.seed',.GlobalEnv)else NULL,seed)
  }
  options(npRmpi.autodispatch.option.sync='invalid');expect_error(check(),"'arg' should be one of",fixed=TRUE)
})

test_that("unverified automatic requests stop before materialization and leave pool usable", {
  if(exists('spawn_mpi_slaves',mode='function')) {
    stopifnot(spawn_mpi_slaves(as.integer(Sys.getenv('N5_WORKERS','1'))))
    on.exit(close_mpi_slaves(force=TRUE),add=TRUE)
  }
  withr::local_options(np.messages=FALSE,npRmpi.autodispatch=TRUE,npRmpi.autodispatch.option.sync='onchange')
  ns<-asNamespace('npRmpi');n<-24L;x<-data.frame(x=seq(.03,.97,length.out=n));z<-data.frame(z=cos(seq_len(n)));y<-sin(4*x$x)+seq_len(n)/100;yd<-data.frame(y=y)
  u<-npudensbw(dat=x,bws=.3,bandwidth.compute=FALSE);d<-npudistbw(dat=x,bws=.3,bandwidth.compute=FALSE)
  rb<-npregbw(xdat=x,ydat=y,bws=.3,bandwidth.compute=FALSE,regtype='ll');cd<-npcdensbw(xdat=x,ydat=yd,bws=c(.3,.3),bandwidth.compute=FALSE);cf<-npcdistbw(xdat=x,ydat=yd,bws=c(.3,.3),bandwidth.compute=FALSE)
  sb<-npscoefbw(xdat=x,ydat=y,zdat=z,bws=.3,bandwidth.compute=FALSE);pb<-npplregbw(xdat=x,ydat=y,zdat=z,bws=matrix(c(.3,.3),2),bandwidth.compute=FALSE);ib<-npindexbw(xdat=cbind(x,z),ydat=y,bws=c(1,.3,.3),method='ichimura',bandwidth.compute=FALSE)
  fit<-npreg(rb,se=FALSE);reference<-fitted(fit);saved<-fit
  # This witness is at the first operation that could force/publish numerical payloads.
  witness<-new.env(parent=emptyenv());witness$n<-0L
  trace('.npRmpi_autodispatch_materialize_call',tracer=function(){witness$n<-witness$n+1L},where=ns,print=FALSE)
  on.exit(untrace('.npRmpi_autodispatch_materialize_call',where=ns),add=TRUE)
  options(npRmpi.autodispatch.option.sync='never')
  calls<-list(quote(npudens(u,tdat=x)),quote(npudist(d,tdat=x)),quote(npreg(rb,txdat=x,tydat=y,se=FALSE)),quote(npcdens(cd,txdat=x,tydat=yd)),quote(npcdist(cf,txdat=x,tydat=yd)),quote(npscoef(sb,txdat=x,tydat=y,tzdat=z)),quote(npplreg(pb,txdat=x,tydat=y,tzdat=z)),quote(npindex(ib,txdat=cbind(x,z),tydat=y)),quote(predict(fit,newdata=x)))
  for(expr in calls){set.seed(7401);seed<-.Random.seed;expect_error(eval(expr),'automatic numerical calls require verified option agreement',fixed=TRUE,info=deparse(expr));expect_identical(.Random.seed,seed)}
  expect_identical(witness$n,0L);expect_identical(fit,saved);expect_identical(fitted(fit),reference)
  # Scalar/small quantile work is local; an explicitly parallel request is not.
  local.q<-npqreg(cf,txdat=x,tydat=yd);expect_true(all(is.finite(fitted(local.q))))
  options(npRmpi.npqreg.parallel.min.eval=1L)
  expect_error(npqreg(cf,txdat=x,tydat=yd),'automatic numerical calls require verified option agreement',fixed=TRUE)
  options(npRmpi.npqreg.parallel.min.eval=NULL)
  expect_identical(witness$n,0L)
  options(npRmpi.autodispatch.option.sync='always');expect_equal(fitted(npreg(rb,se=FALSE)),reference,tolerance=0)
  options(npRmpi.autodispatch.option.sync='onchange');expect_equal(fitted(npreg(rb,se=FALSE)),reference,tolerance=0)
  # An explicit local owner is unaffected by the automatic-call policy.
  options(npRmpi.autodispatch=FALSE,npRmpi.autodispatch.option.sync='never')
  expect_equal(fitted(npreg(rb,se=FALSE)),reference,tolerance=1e-12)
})
