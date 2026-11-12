test_that("private search entry callsites preserve their owner before transport", {
  skip_on_cran()
  result <- npRmpi_run_isolated_contract(c(
    "loadNamespace('npRmpi')",
    "npRmpi::npRmpi.init(nslaves=1L,quiet=TRUE)",
    "on.exit(npRmpi::npRmpi.quit(force=TRUE),add=TRUE)",
    "options(np.messages=FALSE)",
    "ns<-asNamespace('npRmpi'); e<-new.env(parent=ns)",
    "for(n in c('.npregbw_auxiliary','.np_regression_forward','.npregbw_null_impl','.npregbw_default_impl','.npregbw_bandwidth_impl')){f<-get(n,ns);environment(f)<-e;e[[n]]<-f}",
    "e$.npRmpi_autodispatch_call<-function(mc,caller_env,comm=1L,owner.name=NULL) stop(structure(list(message='capture',call=NULL,invocation=mc,owner=owner.name),class=c('owner_capture','error','condition')))",
    "x<-data.frame(x=seq(.05,.95,length.out=24));y<-sin(3*x$x)+cos(1:24)*.1",
    "b<-npRmpi::npregbw(xdat=x,ydat=y,bws=.3,bandwidth.compute=FALSE)",
    "for(id in c('NULL','default','rbandwidth')){ aa<-list(xdat=x,ydat=y);if(id!='NULL')aa$bws<-if(id=='default').3 else b;got<-tryCatch(do.call(e$.npregbw_auxiliary,aa),owner_capture=identity);owner<-paste0('npregbw.',id);stopifnot(inherits(got,'owner_capture'),identical(got$owner,owner),identical(got$invocation[[1]],as.name(owner)),isTRUE(attr(got$invocation,'.npRmpi.bound.call.owner',exact=TRUE)));cat('CALLSITE_CHECKED',id,'\n')}",
    "d<-data.frame(x=x$x,y=y)",
    "for(engine in c('nomad','cell')){set.seed(51);z<-npRmpi::npregbw(y~x,data=d,nomad=(engine=='nomad'),regtype='lp',degree.select='coordinate',search.engine=engine,degree.max=1L,nmulti=1L,nomad.opts=list(MAX_BB_EVAL=8L));stopifnot(inherits(z,'rbandwidth'),is.finite(z$fval))}",
    "cat('SEARCH_CALLSITE_OWNER_PASS\n')"
  ), marker="SEARCH_CALLSITE_OWNER_PASS", timeout=120L)
  if (is.null(result)) skip("installed subprocess library unavailable")
  expect_identical(result$status,0L,info=paste(result$output,collapse="\n"))
  expect_true(result$witnessed)
})
