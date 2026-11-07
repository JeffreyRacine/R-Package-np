# The default suite exercises policy boundaries without compiling custom RNGs.
# Real hidden/exposed generators and public endpoints are qualified in R20.
test_that("signed dependent CDF checks resolved uncached multi-tile replay", {
  ns <- asNamespace("npRmpi")
  infer <- get(".np_udist_joint_inference", ns)
  e <- new.env(parent = environment(infer)); environment(infer) <- e
  e$.np_udist_bc_grid <- function(...) invisible(NULL)
  e$.np_udist_bc_budget <- function(...) list(bytes=0,per.process=0,processes=1)
  e$.np_udist_bc_block_request <- function(x,law,blocklen) list(law=law,length=blocklen)
  e$.np_plot_replay_uniform_state_absent <- function(kind) TRUE
  e$.np_udist_cache_prepare <- function(...) stop("CACHE_OWNER",call.=FALSE)
  e$.np_udist_bc_rows <- function(...) stop("ROWS_OWNER",call.=FALSE)
  state <- list(n=65L,x=data.frame(x=seq_len(65L)),q=1L,s=0L,
    type="fixed",max.bytes=256*1024^2)
  call <- function(m=26L,law="geom",blocklen=1L,calibration="separate-tails",B=97L)
    infer(state,data.frame(x=seq_len(m)),law=law,blocklen=blocklen,
      calibration=calibration,B=B)
  for(law in c("fixed","geom")) for(cal in c("separate-tails","tail-envelope")) {
    expect_error(call(law=law,calibration=cal),"unsafe RNG replay: resolved signed CDF",fixed=TRUE)
    expect_error(call(m=25L,law=law,calibration=cal),"ROWS_OWNER",fixed=TRUE)
    expect_error(call(blocklen=3L,law=law,calibration=cal),"CACHE_OWNER",fixed=TRUE)
  }
  expect_error(call(calibration=NULL),"ROWS_OWNER",fixed=TRUE)
  e$.np_plot_replay_uniform_state_absent <- function(kind) FALSE
  expect_error(call(),"ROWS_OWNER",fixed=TRUE)
})

test_that("signed CDF cache boundaries preserve incumbent selection", {
  plan <- get(".np_udist_cache_plan",asNamespace("npRmpi"))
  f <- function(B=1024L,m=26L,blocklen=3L,max.bytes=Inf)
    plan(2048L,B,m,25L,blocklen,list(bytes=100),max.bytes,626L)
  expect_false(is.null(f(B=1023L)))
  expect_false(is.null(f(B=1024L)))
  expect_null(f(B=1025L))
  expect_null(f(m=25L)); expect_null(f(blocklen=1L))
  need <- f()$bytes
  expect_null(f(max.bytes=need-1))
  expect_false(is.null(f(max.bytes=need)))
  expect_false(is.null(f(max.bytes=need+1)))
})


test_that("unsigned fixed-block CDF refuses only its resolved unsafe tile replay", {
  ns <- asNamespace("npRmpi")
  infer <- get(".np_udist_joint_inference", ns)
  e <- new.env(parent=environment(infer)); environment(infer) <- e
  e$.np_udist_bc_grid <- function(...) invisible(NULL)
  e$.np_udist_bc_budget <- function(...) list(bytes=0,per.process=0,processes=1)
  e$.np_udist_bc_block_request <- function(x,law,blocklen) list(law=law,length=blocklen)
  e$.np_plot_replay_uniform_state_absent <- function(kind) TRUE
  e$.np_udist_bc_rows <- function(...) stop("ROWS_OWNER",call.=FALSE)
  e$.np_udist_gnn_stream <- function(...) stop("GNN_OWNER",call.=FALSE)
  state <- list(n=65L,x=data.frame(x=seq_len(65L)),q=1L,s=0L,
    type="fixed",max.bytes=256*1024^2)
  call <- function(m=26L,law="fixed",blocklen=1L,calibration=NULL)
    infer(state,data.frame(x=seq_len(m)),law=law,blocklen=blocklen,
      calibration=calibration,B=97L)
  for(bl in c(1L,3L)) for(cal in list(NULL,"default","absolute")) {
    for(m in c(26L,33L))
      expect_error(call(m=m,blocklen=bl,calibration=cal),
        "unsafe RNG replay: resolved unsigned fixed-block CDF",fixed=TRUE)
    for(m in c(24L,25L))
      expect_error(call(m=m,blocklen=bl,calibration=cal),"ROWS_OWNER",fixed=TRUE)
  }
  # Numeric GNN has a distinct replay planner; this unit does not reclassify it.
  state$type <- "generalized_nn"
  expect_error(call(),"GNN_OWNER",fixed=TRUE)
  state$q <- 0L
  expect_error(call(),"unsafe RNG replay: resolved unsigned fixed-block CDF",fixed=TRUE)
  state$type <- "fixed"; state$q <- 1L
  e$.np_plot_replay_uniform_state_absent <- function(kind) FALSE
  expect_error(call(),"ROWS_OWNER",fixed=TRUE)
})
