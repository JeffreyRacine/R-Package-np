test_that("raw-state replay refuses a remote custom generator before advance", {
  en <- new.env(parent=environment(.npRmpi_bootstrap_attach_rng_stream))
  f <- .npRmpi_bootstrap_attach_rng_stream
  environment(f) <- en
  en$RNGkind <- function() c("user-supplied", "Inversion", "Rejection")
  en$.npRmpi_has_active_slave_pool <- function() TRUE
  en$.npRmpi_master_only_mode <- function() FALSE
  set.seed(17); before <- .Random.seed; advanced <- FALSE
  tasks <- list(list(start=1L, bsz=3L))
  expect_error(f(tasks, function(task) { advanced <<- TRUE; runif(3) }),
    "cannot reproduce a user-supplied uniform RNG")
  expect_false(advanced)
  expect_identical(.Random.seed, before)
  # Built-in effective streams retain their owner; local seed replay is also uncertifiable.
  en$RNGkind <- function() c("Mersenne-Twister", "Inversion", "Rejection")
  z <- f(tasks, function(task) runif(task$bsz))
  expect_identical(z[[1L]]$rng_state, before)
  expect_identical(attr(z, "rng_final_state"), .Random.seed)
  en$RNGkind <- function() c("user-supplied", "Inversion", "Rejection")
  en$.npRmpi_has_active_slave_pool <- function() FALSE
  en$.npRmpi_master_only_mode <- function() TRUE
  expect_error(f(tasks, function(task) runif(task$bsz)), "user-supplied generator state cannot be certified", fixed=TRUE)
})


test_that("custom normal task replay refuses before querying or advancing workers", {
  en<-new.env(parent=environment(.npRmpi_bootstrap_attach_rng_stream));f<-.npRmpi_bootstrap_attach_rng_stream;environment(f)<-en
  en$RNGkind<-function()c("Mersenne-Twister","user-supplied","Rejection")
  queried<-0L;advanced<-0L;en$.npRmpi_bootstrap_require_normal_rng<-function(){queried<<-queried+1L;invisible(NULL)}
  tasks<-list(list(start=1L,bsz=3L));advance<-function(task){advanced<<-advanced+1L;runif(task$bsz)}
  set.seed(71);before<-.Random.seed
  expect_error(f(tasks,advance,normal.draws=TRUE),"user-supplied generator state cannot be certified",fixed=TRUE)
  expect_identical(queried,0L);expect_identical(advanced,0L);expect_identical(.Random.seed,before)
  uniform<-f(tasks,advance,normal.draws=FALSE)
  expect_null(uniform[[1L]][["normal.kind",exact=TRUE]]);expect_identical(queried,0L);expect_identical(advanced,1L)
})

test_that("custom normal initialization precedes replay state assignment", {
  f <- .npRmpi_bootstrap_use_task_rng
  en <- new.env(parent = environment(f)); environment(f) <- en
  events <- character()
  en$RNGkind <- function(normal.kind) events <<- c(events, normal.kind)
  en$assign <- function(x, value, envir) events <<- c(events, "seed")
  f(list(normal.kind = "user-supplied", rng_state = 1:3))
  expect_identical(events, c("user-supplied", "seed"))
  events <- character(); f(list(rng_state = 1:3))
  expect_identical(events, "seed")
})

test_that("every geometric replay owner declares its indirect normal dependency", {
  ns <- asNamespace("npRmpi")
  owners <- list(
    .np_bias_process_draws=list(F=matrix(1,2,120),B=9L),
    .np_bias_score_student_dispatch=list(F=matrix(1,2,120),active=c(TRUE,TRUE),observed=c(1,1),B=9L,studentize=FALSE),
    .np_udens_b2_boot_dispatch=list(state=list(n=120L),grid=NULL,base=NULL,B=9L),
    .np_udist_bc_block_dispatch=list(F=matrix(1,2,120),raw=NULL,B=9L),
    .np_udist_bc_nn_dispatch=list(state=list(n=120L),grid=NULL,base=NULL,B=9L,corrected=TRUE),
    .np_udist_gnn_draws=list(prepared=list(n=120L),base=1:2,B=9L))
  for (owner in names(owners)) {
    f <- get(owner,ns); e <- new.env(parent=ns); environment(f)<-e
    e$.npRmpi_bootstrap_fanout_enabled<-function(...)TRUE
    e$.npRmpi_bootstrap_chunk_tasks<-function(...)list(list(start=1L,bsz=9L))
    e$.npRmpi_bootstrap_attach_rng_stream<-function(tasks,advance,normal.draws=FALSE) {
      observed <<- normal.draws;stop("captured before advance")
    }
    for (law in c("inid","fixed","geom")) for (bl in c(1L,20L)) {
      observed<-NULL
      expect_error(do.call(f,c(owners[[owner]],list(law=law,blocklen=bl))),"captured before advance")
      expect_identical(observed, law=="geom" && bl>1L, info=paste(owner,law,bl))
    }
  }
  f<-get(".npRmpi_bootstrap_counts_drawer_tasks",ns); e<-new.env(parent=ns);environment(f)<-e
  e$.npRmpi_bootstrap_attach_rng_stream<-function(tasks,advance,normal.draws=FALSE) normal.draws
  drawer<-get(".np_block_counts_drawer",ns)
  for(law in c("fixed","geom"))for(bl in c(1L,19L,20L,21L)) {
    d<-drawer(120L,9L,bl,law)
    expect_identical(f(9L,32L,d),law=="geom" && bl>1L)
  }
  expect_false(f(9L,32L,function(a,b)NULL))
})
