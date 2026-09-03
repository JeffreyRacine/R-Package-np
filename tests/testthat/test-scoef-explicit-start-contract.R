test_that("fixed finalization respects typed zero categorical bandwidths", {
  expect_identical(.npscoef_finalize_bandwidth(c(0,.8),"fixed",24L,
    lower=c(0,1e-8),icon=c(FALSE,TRUE)),c(0,.8))
  expect_identical(.npscoef_finalize_bandwidth(c(.8,0),"fixed",24L,
    lower=c(1e-8,0),icon=c(TRUE,FALSE)),c(.8,0))
  expect_identical(.npscoef_finalize_bandwidth(c(0,.2),"fixed",24L,
    icon=c(FALSE,FALSE)),c(0,.2))
  expect_identical(.npscoef_finalize_bandwidth(c(0,0),"fixed",24L,
    icon=c(FALSE,FALSE)),c(0,0))
  expect_identical(.npscoef_finalize_bandwidth(c(.2,.8),"fixed",24L),c(.2,.8))
  expect_error(.npscoef_finalize_bandwidth(c(0,.8),"fixed",24L),
    "strictly positive",fixed=TRUE)
  expect_error(.npscoef_finalize_bandwidth(c(.2,0),"fixed",24L,icon=c(FALSE,TRUE)),
    "strictly positive",fixed=TRUE)
  expect_error(.npscoef_finalize_bandwidth(c(-.1,.8),"fixed",24L,icon=c(FALSE,TRUE)),
    "categorical bandwidth must be nonnegative",fixed=TRUE)
  expect_error(.npscoef_finalize_bandwidth(c(0,.8),"fixed",24L,icon=FALSE),
    "invalid bandwidth coordinate map",fixed=TRUE)
  expect_error(.npscoef_finalize_bandwidth(c(0,.8),"fixed",24L,icon=c(FALSE,NA)),
    "invalid bandwidth coordinate map",fixed=TRUE)
  expect_error(.npscoef_finalize_bandwidth(c(0,.8),"fixed",24L,
    icon=c(FALSE,TRUE),lower=c(0,1)), "lower bound",fixed=TRUE)
  for(type in c("generalized_nn","adaptive_nn"))
    expect_identical(.npscoef_finalize_bandwidth(c(0,12),type,24L,icon=c(FALSE,TRUE)),
      c(0,12))
})

test_that("explicit zero starts reach the real optimizer and invalid starts do not", {
  expect_true(spawn_mpi_slaves(1L))
  on.exit(close_mpi_slaves(force=TRUE),add=TRUE)
  old <- options(np.messages=FALSE,np.tree=FALSE,np.extendednn=FALSE)
  on.exit(options(old),add=TRUE)
  ns <- asNamespace("npRmpi")
  state <- new.env(parent=emptyenv());state$starts <- list()
  tracer <- substitute({
    optim <- function(par,...) {
      record.state <- STATE
      record.state$starts[[length(record.state$starts)+1L]] <- par
      stats::optim(par,...)
    }
  },list(STATE=state))
  trace("npscoefbw.scbandwidth",where=ns,tracer=tracer,print=FALSE)
  on.exit(untrace("npscoefbw.scbandwidth",where=ns),add=TRUE)
  i <- seq_len(24L)
  x <- data.frame(x1=cos(i/5)+i/24,x2=sin(i/7))
  z1 <- sin(i*sqrt(2))+i/24
  y <- (1+z1/4)*x$x1+(.4+z1^2/10)*x$x2+sin(i*sqrt(5))/20
  z <- data.frame(u=factor(rep(c("a","b"),12L)),z=z1)
  value <- npscoefbw(xdat=x,ydat=y,zdat=z,bws=c(0,12),bwtype="generalized_nn",
    regtype="lc",nomad=FALSE,nmulti=1L,optim.method="Nelder-Mead",
    optim.maxit=120L,optim.maxattempts=1L,random.seed=9303L)
  expect_identical(state$starts[[1L]],c(0,12))
  expect_true(is.finite(value$fval))
  state$starts <- list()
  expect_error(npscoefbw(xdat=x,ydat=y,zdat=z,bws=c(.2,1),bwtype="generalized_nn",
    regtype="lc",nomad=FALSE,nmulti=1L,optim.method="Nelder-Mead"),
    "nearest-neighbor bandwidth",fixed=TRUE)
  expect_length(state$starts,0L)
})
