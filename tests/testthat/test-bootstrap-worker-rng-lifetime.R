test_that("ordinary worker RNG lifetime restores visible state on every unwind", {
  ns <- asNamespace("npRmpi")
  run <- get(".npRmpi_bootstrap_preserve_worker_rng", ns)
  had <- exists(".Random.seed", .GlobalEnv, inherits = FALSE)
  old <- if (had) get(".Random.seed", .GlobalEnv) else NULL
  kinds <- RNGkind()
  on.exit({
    suppressWarnings(do.call(RNGkind, as.list(kinds)))
    if (had) assign(".Random.seed", old, .GlobalEnv) else
      if (exists(".Random.seed", .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir=.GlobalEnv)
  })
  for (kind in c("Mersenne-Twister", "L'Ecuyer-CMRG")) for (present in c(TRUE,FALSE)) {
    RNGkind(kind, "Inversion", "Rejection"); set.seed(314)
    if (!present) rm(".Random.seed", envir=.GlobalEnv)
    entry <- if (present) .Random.seed else NULL
    for (how in c("success", "error", "interrupt", "warning")) {
      f <- function(how) {
        RNGkind("Wichmann-Hill"); set.seed(77); runif(7)
        if (how=="error") stop(structure(list(message="worker failure",call=NULL),class=c("rng_test_error","error","condition")))
        if (how=="interrupt") stop(structure(list(message="worker interrupt",call=NULL),class=c("interrupt","condition")))
        if (how=="warning") warning("worker warning")
        12L
      }
      result <- tryCatch(withCallingHandlers(run(f,how),warning=function(w)invokeRestart("muffleWarning")),
        error=identity,interrupt=identity)
      if (how=="error") expect_s3_class(result,"rng_test_error")
      if (how=="interrupt") expect_s3_class(result,"interrupt")
      if (how%in%c("success","warning")) expect_identical(result,12L)
      expect_identical(RNGkind(),c(kind,"Inversion","Rejection"))
      expect_identical(exists(".Random.seed",.GlobalEnv,inherits=FALSE),present)
      if (present) expect_identical(.Random.seed,entry)
    }
  }
})

test_that("worker lifetime leaves values, body conditions and continuation draws intact", {
  run <- getFromNamespace(".npRmpi_bootstrap_preserve_worker_rng","npRmpi")
  old <- RNGkind(); had <- exists(".Random.seed",.GlobalEnv,inherits=FALSE)
  seed <- if(had) .Random.seed else NULL
  on.exit({suppressWarnings(do.call(RNGkind,as.list(old))); if(had)assign(".Random.seed",seed,.GlobalEnv)
    else if(exists(".Random.seed",.GlobalEnv,inherits=FALSE))rm(".Random.seed",envir=.GlobalEnv)})
  for(kind in c("Mersenne-Twister","L'Ecuyer-CMRG")) {
    RNGkind(kind,"Inversion","Rejection");set.seed(555); expected<-c(runif(5),rnorm(5))
    set.seed(555); value<-run(function(a,b) {set.seed(88);runif(20);list(a=a,b=b)},a=2,b="x")
    expect_identical(value,list(a=2,b="x"));expect_identical(c(runif(5),rnorm(5)),expected)
  }
  RNGkind("Mersenne-Twister","Box-Muller","Rejection");set.seed(555)
  invisible(rnorm(1)); expected<-rnorm(1)
  set.seed(555);invisible(rnorm(1));expect_identical(run(function() 7L),7L)
  expect_identical(rnorm(1),expected)
  RNGkind("Mersenne-Twister","Inversion","Rejection");set.seed(16);entry<-.Random.seed
  caught<-tryCatch(run(function(){set.seed(17);warning("stop from caller");0}),warning=identity)
  expect_s3_class(caught,"warning");expect_identical(.Random.seed,entry)
  suppressWarnings(RNGkind("L'Ecuyer-CMRG","Inversion","Rounding"))
  rm(".Random.seed",envir=.GlobalEnv)
  old.warn <- getOption("warn");options(warn=2)
  on.exit(options(warn=old.warn),add=TRUE)
  expect_identical(run(function(){RNGkind("Mersenne-Twister","Inversion","Rejection");set.seed(2);7L}),7L)
  expect_identical(RNGkind(),c("L'Ecuyer-CMRG","Inversion","Rounding"))
  expect_false(exists(".Random.seed",.GlobalEnv,inherits=FALSE))
})

test_that("preparation separates wrapper bindings and preserves master ownership", {
  prepare <- getFromNamespace(".npRmpi_bootstrap_prepare_worker","npRmpi")
  body <- function(task, extra) list(task=task, extra=extra, worker=worker, comm=comm, preserve=preserve)
  prepared <- prepare(body, list(worker="user worker",comm="user comm",preserve="user preserve"),comm=7L)
  seen <- NULL
  environment(prepared)$mpi.comm.rank <- function(comm) {seen <<- comm; 1L}
  expected <- list(task=3L,extra="value",worker="user worker",comm="user comm",preserve="user preserve")
  expect_identical(prepared(3L,extra="value"),expected)
  expect_identical(seen,7L)
  prepared <- prepare(function(task) {set.seed(task); runif(3)},comm=7L)
  environment(prepared)$mpi.comm.rank <- function(comm) 0L
  set.seed(41);expected <- runif(3);state <- .Random.seed
  set.seed(92);expect_identical(prepared(41),expected)
  expect_identical(.Random.seed,state)
})
