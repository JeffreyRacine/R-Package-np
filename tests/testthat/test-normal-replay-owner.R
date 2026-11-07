test_that("normal replay checks the saved kind and preserves safe decisions", {
  check <- get(".np_bias_normal_replay_check", environment(get(".np_bias_inference", asNamespace("npRmpi"))))
  old <- RNGkind(); on.exit(do.call(RNGkind, as.list(old)), add=TRUE)
  RNGkind("Mersenne-Twister", "Box-Muller"); set.seed(43)
  bm <- .Random.seed
  RNGkind(normal.kind="Inversion"); set.seed(43)
  inv <- .Random.seed; before <- .Random.seed
  for (law in c("gaussian", "geom")) {
    expect_error(check(law, 3L, TRUE, seed=bm), "unsafe RNG replay", fixed=TRUE)
    expect_null(check(law, 3L, FALSE, seed=bm))
    expect_null(check(law, 3L, TRUE, retained=TRUE, seed=bm))
    expect_null(check(law, 3L, TRUE, seed=inv))
  }
  for (law in c("mammen", "rademacher", "inid", "fixed"))
    expect_null(check(law, 3L, TRUE, seed=bm))
  expect_null(check("geom", 1L, TRUE, seed=bm))
  expect_identical(.Random.seed, before)
})

test_that("smoothed reconstruction follows its stream header, not caller kind", {
  ns <- asNamespace("npRmpi"); fn <- function(name) get(name, ns)
  old <- RNGkind(); on.exit(do.call(RNGkind,as.list(old)),add=TRUE)
  x <- data.frame(x=seq(-1,1,length.out=11))
  RNGkind("Mersenne-Twister","Box-Muller");set.seed(42)
  bm <- .Random.seed
  plan <- function(stream,cache=FALSE,layout=NULL,law="inid",blocklen=NULL)
    fn(".np_udens_smoothed_plan")(x,.3,3L,law,blocklen,
      list(count=stream,jitter=stream),2L,cache,layout)
  expect_error(plan(bm),"unsafe RNG replay",fixed=TRUE)
  expect_type(plan(bm,cache=TRUE),"list")
  streams <- fn(".np_udens_calibrated_streams")(bm)
  RNGkind("Mersenne-Twister","Box-Muller");set.seed(42)
  before <- .Random.seed
  expect_type(fn(".np_udens_smoothed_plan")(x,.3,3L,"inid",NULL,
    streams$original,2L,FALSE),"list")
  expect_identical(.Random.seed,before)
})

# Resolve the master-owned decision without allocating an MPI pool in a unit test.
# Actual worker/cache execution and transport refusals have separate MPI gates.
test_that("CDF master refuses before an uncached normal producer", {
  ns <- asNamespace("npRmpi");infer <- get(".np_udist_joint_inference",ns)
  e <- new.env(parent=environment(infer));environment(infer) <- e
  e$.np_udist_bc_grid <- function(...) invisible(NULL)
  e$.np_udist_bc_budget <- function(...) list(bytes=0,per.process=0,processes=1)
  e$.np_udist_bc_block_request <- function(x,law,blocklen) list(law=law,length=blocklen)
  e$.np_udist_bc_rows <- function(...) stop("ROWS_OWNER",call.=FALSE)
  old <- RNGkind();on.exit(do.call(RNGkind,as.list(old)),add=TRUE)
  state <- list(n=65L,x=data.frame(x=seq_len(65L)),q=1L,s=0L,
    type="fixed",max.bytes=256*1024^2)
  call <- function(m=26,law="geom",blocklen=3L,bands="simultaneous") {
    RNGkind("Mersenne-Twister","Box-Muller");set.seed(43)
    infer(state,data.frame(x=seq_len(m)),law=law,blocklen=blocklen,B=49L,bands=bands)
  }
  expect_error(call(),"unsafe RNG replay",fixed=TRUE)
  expect_error(call(law="gaussian",blocklen=NULL),"unsafe RNG replay",fixed=TRUE)
  expect_error(call(m=25),"ROWS_OWNER",fixed=TRUE)
  expect_error(call(blocklen=1L),"ROWS_OWNER",fixed=TRUE)
  expect_error(call(law="gaussian",blocklen=NULL,bands="pointwise"),"ROWS_OWNER",fixed=TRUE)
})

test_that("regression checks resolved block and family replay before row work", {
  ns <- asNamespace("npRmpi");infer <- get(".np_bias_inference",ns)
  e <- new.env(parent=environment(infer));environment(infer) <- e
  e$.np_bias_regression_rows <- function(...) stop("ROWS_OWNER",call.=FALSE)
  state <- list(x=data.frame(x=seq_len(41)),residual_raw=sin(seq_len(41)))
  old <- RNGkind();on.exit(do.call(RNGkind,as.list(old)),add=TRUE)
  call <- function(m=26,law="geom",L=3L,targets="mean",calibration="absolute") {
    RNGkind("Mersenne-Twister","Box-Muller");set.seed(43)
    infer(state,data.frame(x=seq_len(m)),targets=targets,law=law,
      blocklen=if(law=="geom")L else NULL,B=49L,calibration=calibration)
  }
  expect_error(call(),"unsafe RNG replay",fixed=TRUE)
  expect_error(call(targets=c("mean","d1")),"unsafe RNG replay",fixed=TRUE)
  expect_error(call(calibration="separate-tails"),"unsafe RNG replay",fixed=TRUE)
  expect_error(call(m=25),"ROWS_OWNER",fixed=TRUE)
  expect_error(call(L=1L),"ROWS_OWNER",fixed=TRUE)
  expect_error(call(law="gaussian"),"ROWS_OWNER",fixed=TRUE)
  # Stand-in for the resolved uncached representation; real custom generators
  # and public endpoints are exercised in the fresh-process forensic gate.
  e$.np_bias_gaussian_cache <- function(...) NULL
  expect_error(call(law="gaussian"),"unsafe RNG replay",fixed=TRUE)
  expect_error(call(m=25,law="gaussian"),"ROWS_OWNER",fixed=TRUE)
})
