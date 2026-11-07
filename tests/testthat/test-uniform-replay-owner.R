test_that("resolved replay refuses known missing uniform payload without draws", {
  ns <- asNamespace("npRmpi")
  check <- get(".np_bias_replay_check", environment(get(".np_bias_inference", ns)))
  old <- RNGkind(); had <- exists(".Random.seed", .GlobalEnv, inherits=FALSE)
  before <- get0(".Random.seed", .GlobalEnv, inherits=FALSE)
  on.exit({do.call(RNGkind,as.list(old));if(had)assign(".Random.seed",before,.GlobalEnv) else if(exists(".Random.seed",.GlobalEnv,inherits=FALSE))rm(".Random.seed",envir=.GlobalEnv)},add=TRUE)
  RNGkind("Mersenne-Twister","Inversion");set.seed(42); seed <- .Random.seed
  # Header fields are established R RNG identifiers; actual generators have a
  # separate fresh-process public proof. No fake seed is installed here.
  for(law in c("gaussian","mammen","rademacher","inid","fixed","geom")) {
    expect_error(check(law,3L,TRUE,seed=10405L),"custom uniform state absent",fixed=TRUE)
    expect_null(check(law,3L,FALSE,seed=10405L))
    expect_null(check(law,3L,TRUE,retained=TRUE,seed=10405L))
    expect_null(check(law,3L,TRUE,seed=c(10405L,17L)))
    expect_null(check(law,3L,TRUE,seed=seed))
  }
  # Normal callbacks do not necessarily use uniform state; that opaque class
  # stays separate, while count draws still require it.
  expect_null(check("gaussian",NULL,TRUE,seed=10305L))
  expect_error(check("inid",NULL,TRUE,seed=10305L),"custom uniform state absent",fixed=TRUE)
  expect_error(check("gaussian",NULL,TRUE,seed=10203L),"Box-Muller",fixed=TRUE)
  expect_identical(.Random.seed,seed)
})

test_that("density inid checks replay before contribution work", {
  ns <- asNamespace("npRmpi"); f <- get(".np_udens_bc_inid",ns)
  e<-new.env(parent=environment(f));environment(f)<-e
  e$.np_udens_bc_contribution<-function(...)stop("CONTRIBUTION_OWNER",call.=FALSE)
  grid<-data.frame(x=seq_len(26))
  expect_error(f(list(),grid,.05,49L,25L,TRUE,10405L),"custom uniform state absent",fixed=TRUE)
  expect_error(f(list(),grid[1:25,,drop=FALSE],.05,49L,25L,TRUE,10405L),"CONTRIBUTION_OWNER",fixed=TRUE)
})

test_that("GNN replay uses actual marginal-profile blocks", {
  ns<-asNamespace("npRmpi");f<-get(".np_udist_gnn_stream",ns)
  e<-new.env(parent=environment(f));environment(f)<-e
  e$.np_udist_gnn_budget<-function(...)list(marginal.tile=1L)
  e$.np_udist_gnn_main<-function(...)list(center=0,original.center=0)
  e$.np_udist_gnn_blocks<-function(...)list(list(),list())
  e$.np_udist_gnn_block<-function(...)stop("BLOCK_OWNER",call.=FALSE)
  expect_error(f(list(),data.frame(x=0),"inid",.05,49L,TRUE,NULL,10405L),"custom uniform state absent",fixed=TRUE)
  e$.np_udist_gnn_blocks<-function(...)list(list())
  expect_error(f(list(),data.frame(x=0),"inid",.05,49L,TRUE,NULL,10405L),"BLOCK_OWNER",fixed=TRUE)
})

test_that("smoothed categorical jitter checks its own saved state", {
  ns<-asNamespace("npRmpi");f<-get(".np_udens_smoothed_plan",ns)
  e<-new.env(parent=environment(f));environment(f)<-e
  e$.np_udens_smoothed_samples<-function(...)stop("SAMPLES_OWNER",call.=FALSE)
  streams<-list(count=c(10407L,rep(1L,6)),jitter=10405L)
  layout<-list(numeric=FALSE,category=list(list(index=1L,ordered=FALSE,lambda=0)))
  call<-function(cache=FALSE)f(data.frame(x=1:3),.2,3L,"inid",NULL,streams,2L,cache,layout)
  expect_error(call(),"custom uniform state absent",fixed=TRUE)
  expect_error(call(TRUE),"SAMPLES_OWNER",fixed=TRUE)
  streams$jitter<-c(10407L,rep(1L,6))
  expect_error(call(),"SAMPLES_OWNER",fixed=TRUE)
})
