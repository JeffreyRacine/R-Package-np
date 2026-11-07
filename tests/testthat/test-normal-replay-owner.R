test_that("normal replay checks the saved kind and preserves safe decisions", {
  check <- get(".np_bias_normal_replay_check", environment(get(".np_bias_inference", asNamespace("np"))))
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

test_that("CDF owners refuse uncached normal replay and retain caches", {
  ns <- asNamespace("np"); fn <- function(name) get(name, ns)
  old <- RNGkind(); on.exit(do.call(RNGkind,as.list(old)),add=TRUE)
  x <- data.frame(x=sin(seq_len(41)))
  bw <- npudistbw(dat=x,bws=.4,bandwidth.compute=FALSE)
  state <- fn(".np_udist_bc_prepare")(bw,x)
  grid <- data.frame(x=seq(-.8,.8,length.out=26))
  reset <- function() {RNGkind("Mersenne-Twister","Box-Muller");set.seed(83)}
  call <- function(m=26, law="geom", blocklen=3L, calibration=NULL, bands="simultaneous") {
    reset()
    fn(".np_udist_joint_inference")(state,grid[seq_len(m),,drop=FALSE],law=law,
      blocklen=blocklen,B=49L,calibration=calibration,bands=bands)
  }
  expect_error(call(),"unsafe RNG replay",fixed=TRUE)
  expect_error(call(law="gaussian",blocklen=NULL),"unsafe RNG replay",fixed=TRUE)
  expect_type(call(m=25),"list")
  expect_type(call(blocklen=1L),"list")
  expect_type(call(law="gaussian",blocklen=NULL,bands="pointwise"),"list")
  expect_type(call(calibration="separate-tails"),"list")
  reset()
  expect_error(fn(".np_udist_bc_inference")(state,grid,law="geom",blocklen=3L,B=49L),
    "unsafe RNG replay",fixed=TRUE)
  reset()
  expect_type(fn(".np_udist_bc_inference")(state,grid,law="gaussian",B=49L),"list")
  larger <- data.frame(x=seq(-.8,.8,length.out=33))
  reset()
  expect_error(fn(".np_udist_bc_inference")(state,larger,law="gaussian",B=49L),
    "unsafe RNG replay",fixed=TRUE)
  # The GNN planner tiles marginal coordinates and categorical profiles separately.
  bw <- npudistbw(dat=x,bws=12L,bwtype="generalized_nn",bandwidth.compute=FALSE)
  state <- fn(".np_udist_bc_prepare")(bw,x)
  reset()
  expect_error(fn(".np_udist_gnn_stream")(state,grid,"geom",.05,49L,TRUE,3L,.Random.seed),
    "unsafe RNG replay",fixed=TRUE)
  reset()
  expect_type(fn(".np_udist_gnn_stream")(state,grid[1:2,,drop=FALSE],"geom",.05,49L,TRUE,3L,.Random.seed),"list")
})

test_that("density owners use factor, family and retained-draw boundaries", {
  ns <- asNamespace("np"); fn <- function(name) get(name, ns)
  old <- RNGkind(); on.exit(do.call(RNGkind,as.list(old)),add=TRUE)
  x <- data.frame(x=sin(seq_len(41)))
  bw <- npudensbw(dat=x,bws=.4,bandwidth.compute=FALSE)
  state <- fn(".np_udens_bc_prepare")(bw,x)
  grid <- data.frame(x=seq(-.8,.8,length.out=26))
  reset <- function() {RNGkind("Mersenne-Twister","Box-Muller");set.seed(83)}
  reset()
  expect_error(fn(".np_udens_bc_factor")(state,grid,.05,49L,25L,TRUE,.Random.seed),
    "unsafe RNG replay",fixed=TRUE)
  reset()
  expect_type(fn(".np_udens_bc_factor")(state,grid[1:25,,drop=FALSE],.05,49L,25L,TRUE,.Random.seed),"list")
  reset()
  expect_error(fn(".np_udens_b2_bootstrap")(state,grid,law="geom",blocklen=3L,B=49L),
    "unsafe RNG replay",fixed=TRUE)
  reset()
  expect_type(fn(".np_udens_b2_bootstrap")(state,grid[1:2,,drop=FALSE],law="geom",blocklen=3L,B=49L),"list")
  reset()
  expect_error(fn(".np_udens_reference_inference")(state,grid,law="geom",blocklen=3L,B=49L),
    "unsafe RNG replay",fixed=TRUE)
  reset()
  expect_type(fn(".np_udens_reference_inference")(state,grid[1:2,,drop=FALSE],law="geom",blocklen=3L,B=49L),"list")
  # A second scale pass replays even when the whole grid fits in one tile.
  reset()
  expect_error(fn(".np_udens_scale_inference")(state,grid[1:2,,drop=FALSE],law="geom",blocklen=3L,B=49L),
    "unsafe RNG replay",fixed=TRUE)
  reset()
  expect_type(fn(".np_udens_scale_inference")(state,grid,law="geom",blocklen=3L,B=49L),"list")
  d8 <- fn(".np_udens_d8_prepare")(bw,x)
  reset()
  expect_error(fn(".np_udens_d8_inference")(d8,grid,B=49L),"unsafe RNG replay",fixed=TRUE)
  reset()
  expect_type(fn(".np_udens_d8_inference")(d8,grid[1:2,,drop=FALSE],B=49L),"list")
})

test_that("smoothed reconstruction follows its stream header, not caller kind", {
  ns <- asNamespace("np"); fn <- function(name) get(name, ns)
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

test_that("regression checks resolved block and family replay before row work", {
  ns <- asNamespace("np");infer <- get(".np_bias_inference",ns)
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
