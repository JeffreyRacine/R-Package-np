test_that("significance formula handoff consumes a prepared NA policy once", {
  pkg <- environmentName(environment(npsigtest))
  if (pkg == "npRmpi") {
    if (!spawn_mpi_slaves()) skip("MPI slaves unavailable")
    on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  }
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  i <- seq_len(30L)
  d <- data.frame(x1 = sin(i/2), x2 = cos(i/3))
  d$y <- d$x1^2 + .5*d$x2 + sin(5*i)/4
  d$x1[4L] <- NA; d$y[9L] <- NA
  keep <- which(!is.na(d$x1) & d$x1 > -.5 & complete.cases(d))
  reference <- d[keep, , drop = FALSE]
  values <- function(x) x[c("In", "In.bootstrap", "P", "ixvar", "boot.num",
    "pivot.effective", "bootstrap.executed", "rejectNum")]
  once <- function(...) npsigtest(...)
  twice <- function(...) once(...)
  for (route in c("positional", "named", "formula", "one", "two")) {
    factory.hits <- policy.hits <- 0L
    policy <- function(x) { policy.hits <<- policy.hits + 1L; stats::na.omit(x) }
    factory <- function() { factory.hits <<- factory.hits + 1L; policy }
    set.seed(8301)
    got <- switch(route,
      positional = npsigtest(y~x1+x2, data=d, subset=x1 > -.5,
        na.action=factory(), nmulti=1L, B=9L, pivot=FALSE),
      named = npsigtest(bws=y~x1+x2, data=d, subset=x1 > -.5,
        na.action=factory(), nmulti=1L, B=9L, pivot=FALSE),
      formula = npsigtest(formula=y~x1+x2, data=d, subset=x1 > -.5,
        na.action=factory(), nmulti=1L, B=9L, pivot=FALSE),
      one = once(y~x1+x2, data=d, subset=x1 > -.5,
        na.action=factory(), nmulti=1L, B=9L, pivot=FALSE),
      two = twice(y~x1+x2, data=d, subset=x1 > -.5,
        na.action=factory(), nmulti=1L, B=9L, pivot=FALSE))
    seed <- .Random.seed
    expect_identical(factory.hits, 1L, info=route)
    expect_identical(policy.hits, 1L, info=route)
    set.seed(8301)
    expected <- npsigtest(xdat=reference[c("x1","x2")], ydat=reference$y,
                         nmulti=1L, B=9L, pivot=FALSE)
    expect_identical(.Random.seed, seed, info=route)
    expect_equal(values(got), values(expected), tolerance=0, info=route)
    expect_identical(got$bws$bw, expected$bws$bw, info=route)
    expect_false(grepl("frame.state|<environment", paste(deparse(got$call),collapse=" ")))
  }
  factory.hits <- policy.hits <- 0L
  set.seed(8305)
  got <- npsigtest(y~x1+x2,data=d,subset=x1 > -.5,na.action=factory(),B=9L)
  seed <- .Random.seed
  expect_identical(factory.hits,1L)
  expect_identical(policy.hits,1L)
  set.seed(8305)
  expected <- npsigtest(xdat=reference[c("x1","x2")],ydat=reference$y,B=9L)
  expect_identical(.Random.seed,seed)
  expect_equal(values(got),values(expected),tolerance=0)
  expect_identical(got$bws$bw,expected$bws$bw)
})

test_that("significance replacement data owns its NA policy and preserves refusals", {
  pkg <- environmentName(environment(npsigtest))
  if (pkg == "npRmpi") {
    if (!spawn_mpi_slaves()) skip("MPI slaves unavailable")
    on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  }
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  i <- seq_len(30L)
  d <- data.frame(x1=sin(i/2),x2=cos(i/3))
  d$y <- d$x1^2+.5*d$x2+sin(5*i)/4
  bw <- npregbw(y~x1+x2,data=d,bws=c(.5,.5),bandwidth.compute=FALSE)
  before <- bw
  dna <- d; dna$x1[4L] <- NA; dna$y[9L] <- NA
  reference <- dna[complete.cases(dna),,drop=FALSE]
  values <- function(x) x[c("In", "In.bootstrap", "P", "ixvar", "boot.num",
    "pivot.effective", "bootstrap.executed", "rejectNum")]
  for (method in c("iid","wild","wild-rademacher","pairwise","typeII")) {
    factory.hits <- policy.hits <- 0L
    policy <- function(x) {policy.hits <<- policy.hits+1L; stats::na.omit(x)}
    factory <- function() {factory.hits <<- factory.hits+1L; policy}
    boot.method <- if(method=="typeII") "iid" else method
    boot.type <- if(method=="typeII") "II" else "I"
    set.seed(8302)
    got <- npsigtest(bw,data=dna,na.action=factory(),B=9L,
                     boot.method=boot.method,boot.type=boot.type)
    seed <- .Random.seed
    expect_identical(factory.hits,1L,info=method)
    expect_identical(policy.hits,1L,info=method)
    set.seed(8302)
    expected <- npsigtest(bw,xdat=reference[c("x1","x2")],ydat=reference$y,
                          B=9L,boot.method=boot.method,boot.type=boot.type)
    expect_identical(.Random.seed,seed,info=method)
    expect_equal(values(got),values(expected),tolerance=0,info=method)
    expect_identical(bw,before)
  }
  for (policy in list(stats::na.omit,"na.omit",NULL)) {
    set.seed(8303)
    got <- npsigtest(bw,data=d,na.action=policy,B=9L)
    seed <- .Random.seed
    set.seed(8303)
    expected <- npsigtest(bw,B=9L)
    expect_identical(.Random.seed,seed)
    expect_equal(values(got),values(expected),tolerance=0)
  }
  expect_error(npsigtest(bw,na.action=stats::na.omit,B=9L),
               "na.action override is not supported",fixed=TRUE)
  expect_error(npsigtest(bw,data=NULL,na.action=stats::na.omit,B=9L),
               "na.action override is not supported",fixed=TRUE)
  expect_error(npsigtest(bw,xdat=d[c("x1","x2")],ydat=d$y,
                         na.action=stats::na.omit,B=9L),
               "na.action override is not supported",fixed=TRUE)
  set.seed(8304)
  seed <- .Random.seed
  expect_error(npsigtest(bw,data=dna,na.action=stats::na.fail,B=9L),
               "missing values",fixed=TRUE)
  expect_identical(.Random.seed,seed)
  expect_identical(bw,before)
})

test_that("replacement without an NA override retains lazy option forwarding", {
  pkg <- environmentName(environment(npsigtest))
  method <- getS3method("npsigtest", "formula")
  d <- data.frame(y=seq_len(5L),x=seq_len(5L))
  mf <- model.frame(y~x,data=d)
  b <- list(call=quote(npregbw(y~x,data=d)),formula=y~x,terms=terms(y~x))
  class(b) <- "rbandwidth"
  local_mocked_bindings(
    .np_bws_formula_model_frame = function(...) mf,
    npsigtest = function(xdat, ydat, bws, ..., B=399) {
      values <- list(...)
      list(options=values,B=B)
    }, .package=pkg)
  hits <- 0L
  inert <- quote({hits <- hits+1L; 1})
  result <- method(b,data=d,opt=inert,B=9L)
  expect_identical(hits,0L)
  expect_identical(result$options$opt,inert)
  expect_error(method(b,data=d,na.action=stats::na.omit,opt=inert,B=9L),
               "executable option value",fixed=TRUE)
  expect_identical(hits,0L)
})
