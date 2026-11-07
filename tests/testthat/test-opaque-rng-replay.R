test_that("opaque custom replay refuses only participating uncached streams", {
  ns<-asNamespace("npRmpi");f<-get(".np_bias_replay_check",ns)
  set.seed(8001);before<-.Random.seed;count<-0L
  for(uniform in 0:7)for(normal in 0:5)for(law in c("gaussian","inid","fixed","geom","mammen","rademacher"))for(bl in c(1L,4L)) {
    seed<-as.integer(c(uniform+100L*normal+10000L,rep(123,7)))
    consumes.normal<-law=="gaussian" || (law=="geom" && bl>1L)
    unsafe<-uniform==5L || (consumes.normal && normal%in%c(2L,3L))
    if(unsafe)expect_error(f(law,bl,TRUE,seed=seed),"unsafe RNG replay",fixed=TRUE) else expect_null(f(law,bl,TRUE,seed=seed))
    expect_null(f(law,bl,FALSE,seed=seed));expect_null(f(law,bl,TRUE,retained=TRUE,seed=seed));count<-count+1L
  }
  expect_identical(count,576L);expect_identical(.Random.seed,before)
  # Previously diagnosed absent-uniform and Box-Muller errors retain precedence.
  expect_error(f("inid",1L,TRUE,seed=10405L),"custom uniform state absent",fixed=TRUE)
  expect_error(f("gaussian",1L,TRUE,seed=as.integer(c(10205,1:3))),"Box-Muller normal state",fixed=TRUE)
})


test_that("manual replay retains its common refusal before RNG inspection", {
  ns<-asNamespace("npRmpi");env<-new.env(parent=ns)
  env$.npRmpi_autodispatch_called_from_bcast<-function()TRUE
  fanout<-get(".npRmpi_bootstrap_fanout_enabled",ns);environment(fanout)<-env
  env$.npRmpi_bootstrap_fanout_enabled<-fanout
  guard<-get(".npRmpi_replay_context_preflight",environment(get(".np_bias_replay_check",ns)));environment(guard)<-env
  env$.npRmpi_replay_context_preflight<-guard
  env$RNGkind<-function(...)stop("RNG inspection must not occur")
  message<-"MPI corrected process cannot run inside mpi.bcast.cmd context"
  expect_error(guard(),message,fixed=TRUE)
  replay<-get(".np_bias_replay_check",ns);environment(replay)<-env
  expect_error(replay("gaussian",NULL,TRUE),message,fixed=TRUE)
  expect_silent(replay("gaussian",NULL,FALSE))
  expect_silent(replay("gaussian",NULL,TRUE,retained=TRUE))
  attach<-get(".npRmpi_bootstrap_attach_rng_stream",ns);environment(attach)<-env
  expect_identical(attach(list(),NULL),list())
  expect_error(attach(list(list(start=1L,bsz=2L)),function(...)stop("advance must not occur")),message,fixed=TRUE)
  factor<-get(".np_plot_replay_factor_preflight",ns);environment(factor)<-env
  b<-structure(list(type="fixed",ncon=1L,nuno=0L,nord=0L,ckertype="gaussian",ckerorder=2L,ckerbound="none"),class="dbandwidth")
  d<-list(plot.errors.method="asymptotic",plot.errors.type="simultaneous",neval=33L,bias.control=np_udist_bias_control("cdf"))
  expect_error(factor(get(".np_plot_dbandwidth_engine",ns),b,d),message,fixed=TRUE)
  d$neval<-32L;expect_silent(factor(get(".np_plot_dbandwidth_engine",ns),b,d))
})
