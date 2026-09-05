.test_cdf_raw_result_contract <- function(pkg,degree.search,refine) {
  ns <- asNamespace(pkg)
  old <- options(np.messages=FALSE,np.extendednn=FALSE,np.tree=FALSE)
  on.exit(options(old),add=TRUE)
  n <- 24L
  x <- data.frame(x=sin(seq_len(n)*sqrt(2))+seq_len(n)/n)
  y <- sin(seq_len(n)*sqrt(3))+seq_len(n)/(2*n)
  points <- list(c(8,8),c(18,18))
  degrees <- if(degree.search && !refine) c(0L,1L) else c(0L,0L)
  raw.fn <- get(".npcdistbw_eval_only",ns)
  make <- function(p,d) npcdistbw(xdat=x,ydat=y,bwtype="adaptive_nn",
    bws=p,regtype="lp",degree=d,nomad=FALSE,bandwidth.compute=FALSE)
  raw <- vapply(1:2,function(i) raw.fn(xdat=x,ydat=y,
    bws=make(points[[i]],degrees[[i]]),ngrid=7L,
    invalid.penalty="dbmax")$objective,numeric(1))
  expect_true(all(is.finite(raw) & raw < .Machine$double.xmax))
  expect_false(identical(raw[[1L]],raw[[2L]]))
  winner <- which.min(raw)
  expected.point <- points[[winner]]
  expected.raw <- raw[[winner]]
  expected.degree <- degrees[[winner]]
  fake <- -raw-10
  if(refine) {
    loser <- which.max(raw)
    points <- rep(points[loser],2L)
    fake <- rep(raw[[loser]],2L)
  }
  native.calls <- raw.calls <- hot.calls <- 0L
  native.mock <- function(prep,x0,...) {
    native.calls <<- native.calls+1L
    i <- native.calls
    p <- c(points[[i]],if(degree.search) degrees[[i]])
    list(status=0L,result_status=0L,nomad_run_flag=0L,
      objective=fake[[i]],official_objective=fake[[i]],solution=p,best_point=p,
      best_degree=if(degree.search) degrees[[i]] else integer(),
      first_degree=if(degree.search) degrees[[i]] else integer(),
      first_objective=fake[[i]],best_fval=fake[[i]],best_num.feval=1,
      best_num.feval.fast=0,best_num.feval.guarded=0,best_num.feval.invalid=0,
      total_num.feval=1,total_num.feval.fast=0,total_num.feval.guarded=0,
      total_num.feval.invalid=0,compiled_callback_calls=1L,
      compiled_callback_failures=0L,crs_callback_evaluations=1L,
      blackbox_evaluations=1L,cache_hits=0L,cache_size=0L,
      total_evaluations=1L,iterations=1L,message="controlled native result")
  }
  raw.mock <- function(...) { raw.calls <<- raw.calls+1L; raw.fn(...) }
  hot.mock <- function(...) {
    hot.calls <<- hot.calls+1L
    b <- make(expected.point,0L)
    b$fval <- 100
    b$ifval <- 777
    b$fval.history <- c(777,100)
    b$num.feval <- 1
    b$num.feval.fast <- 0
    b
  }
  bindings <- list(npNomadNativeSearchConditionalDistribution=native.mock,
    .npcdistbw_eval_only=raw.mock,.npcdistbw_run_fixed_degree=hot.mock,
    .package=pkg,.env=environment())
  if(pkg=="npRmpi")
    bindings$.npcdistbw_run_fixed_degree_collective <- hot.mock
  do.call(testthat::local_mocked_bindings,bindings)
  controls <- if(degree.search) list(nomad=TRUE,
    search.engine=if(refine) "nomad+powell" else "nomad",
    degree.min=0L,degree.max=1L,degree.start=0L,degree.verify=FALSE)
    else list(nomad=FALSE,degree=0L,bwsolver=if(refine) "mads+powell" else "mads")
  b <- do.call(npcdistbw,c(list(xdat=x,ydat=y,bwtype="adaptive_nn",
    regtype="lp",nmulti=if(refine) 1L else 2L,ngrid=7L,random.seed=42),controls))
  expect_identical(as.numeric(c(b$ybw,b$xbw)),expected.point)
  expect_identical(b$fval,expected.raw)
  expect_equal(b$degree,expected.degree)
  expect_identical(native.calls,if(refine) 1L else 2L)
  expect_identical(raw.calls,3L)
  expect_identical(hot.calls,if(refine) 1L else 0L)
  expect_identical(b$ifval,if(refine) 777 else expected.raw)
  expect_identical(b$fval.history,if(refine) c(777,100) else expected.raw)
}

.test_cdf_exact_powell_endpoint <- function(pkg) {
  old <- options(np.messages=FALSE,np.extendednn=FALSE,np.tree=FALSE)
  on.exit(options(old),add=TRUE)
  i <- 1:48
  x <- data.frame(x=sin(i*sqrt(2))+i/48,z=factor(rep(c("a","b"),24)))
  y <- x$x+.3*(x$z=="b")+sin(i*sqrt(3))/2
  b <- npcdistbw(xdat=x,ydat=y,bwtype="fixed",regtype="lp",nomad=TRUE,
    search.engine="nomad+powell",degree.min=0L,degree.max=1L,
    degree.start=1L,degree.verify=FALSE,nmulti=1L,itmax=12L,
    powell.remin=FALSE,nomad.opts=list(MAX_BB_EVAL=24L),
    random.seed=42L,ngrid=7L)
  raw <- get(".npcdistbw_eval_only",asNamespace(pkg))(
    xdat=x,ydat=y,bws=b,ngrid=7L,invalid.penalty="dbmax")$objective
  expect_true(is.finite(raw) && abs(raw)<.Machine$double.xmax)
  expect_identical(b$fval,raw)
}

test_that("local CDF Powell certification retains the exact bandwidth point", {
  .npRmpi_with_local_regression(.test_cdf_exact_powell_endpoint("npRmpi"))
})

test_that("local native CDF selection publishes already certified raw objectives", {
  for(degree.search in c(FALSE,TRUE)) for(refine in c(FALSE,TRUE))
    .npRmpi_with_local_regression(.test_cdf_raw_result_contract("npRmpi",degree.search,refine))
})
