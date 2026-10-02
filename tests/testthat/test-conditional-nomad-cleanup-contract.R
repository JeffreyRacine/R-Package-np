test_that("conditional NOMAD terminal errors preserve same-process recovery", {
  skip_if_not_installed("crs")
  pkg <- getNamespaceName(environment(npcdensbw))
  old <- options(np.messages=FALSE, np.tree=FALSE)
  on.exit(options(old), add=TRUE)
  loadNamespace("crs")
  quadratic <- function() .Call("crs_nomad_native_test_solve", list(
    x0=c(.1,.2), lower=c(-1,-1), upper=c(1,1), input_type=c(0L,0L),
    output_type=0L, max_eval=12L, random_seed=42L), PACKAGE="crs")
  expected <- quadratic()
  expect_identical(expected$status, 0L)
  recover <- function() {
    got <- quadratic()
    expect_identical(got$status, 0L)
    expect_identical(got$solution, expected$solution)
    expect_identical(got$objective, expected$objective)
  }
  # Retain the original small-n shortcut ambiguity, including its public text.
  for (visible in c(FALSE, TRUE)) {
    options(np.messages=visible)
    set.seed(20260322)
    dat <- data.frame(x=sort(runif(14)), y=sort(runif(14)))
    expect_error(npcdensbw(y~x, data=dat, nomad=TRUE, degree.max=1L, nmulti=1L),
      "conditional bandwidth search stopped: ambiguous numerical rank; row 14; bandwidth/scale factors")
    recover()
  }
  # Original CDF progress fixture; successful progress coverage is separate.
  options(np.messages=FALSE)
  set.seed(20260401)
  x <- sort(runif(18)); z <- sort(runif(18))
  y <- sin(2*pi*x)+rnorm(18,sd=.05)
  expect_error(npcdistbw(y~x, data=data.frame(x,y), regtype="lp",
    degree.select="coordinate", search.engine="nomad+powell",
    degree.min=0L, degree.max=1L, degree.verify=FALSE, bwtype="fixed",
    bwmethod="cv.ls", nmulti=2L, ngrid=30L),
    "conditional bandwidth search stopped: ambiguous numerical rank; row 18; bandwidth/scale factors")
  recover()
  # Actual conditional native work still succeeds after both failure owners.
  set.seed(952)
  dat <- data.frame(x=seq(-1,1,length.out=48),y=rnorm(48))
  for (fun in list(npcdensbw,npcdistbw)) {
    bw <- fun(y~x, data=dat, regtype="lp", degree=0L,
      bwsolver="mads", bws=c(.5,.8), nmulti=1L,
      nomad.opts=list(MAX_BB_EVAL=12L))
    expect_true(all(is.finite(c(bw$xbw,bw$ybw))))
  }
  recover()
})


test_that("conditional density preparation errors do not poison later penalties", {
  old <- options(np.messages=FALSE, np.tree=FALSE)
  on.exit(options(old), add=TRUE)
  ns <- asNamespace(getNamespaceName(environment(npcdensbw)))
  prepare <- get("npPreparedObjectivePrepareConditionalDensity", ns)
  destroy <- get("npPreparedObjectiveDestroyConditionalDensity", ns)
  evaluate <- get("npPreparedObjectiveEvalConditionalDensity", ns)
  args_for <- function(dat, degree, bandwidth) {
    bw <- npcdensbw(y~x, data=dat, regtype="lp", degree=degree,
                   bws=bandwidth, bandwidth.compute=FALSE)
    args <- get(".npcdensbw_prepared_prepare_args", ns)(
      dat["x"], dat["y"], bw, start.bw=bandwidth, invalid.penalty="baseline")
    names(args)[names(args)=="penalty_mode"] <- "penalty.mode"
    names(args)[names(args)=="penalty_multiplier"] <- "penalty.multiplier"
    args
  }
  set.seed(952)
  healthy <- data.frame(x=seq(-1,1,length.out=48), y=rnorm(48))
  healthy_args <- args_for(healthy, 0L, c(.5,.8))
  penalty <- function() {
    on.exit(destroy(), add=TRUE)
    expect_true(as.logical(do.call(prepare, healthy_args)))
    evaluate(c(-1,-1), 0L)
  }
  before <- penalty()
  set.seed(20260322)
  bad <- data.frame(x=sort(runif(14)), y=sort(runif(14)))
  bandwidth <- c(.048130207494902745,.025634352719319447)
  bad_args <- args_for(bad, 1L, bandwidth)
  expect_error(do.call(prepare, bad_args),
    "conditional bandwidth search stopped: ambiguous numerical rank; row 14;")
  expect_identical(penalty(), before)
  expect_error(npcdensbw(y~x, data=bad, regtype="lp", degree=1L,
    bws=bandwidth, bwsolver="mads", nmulti=1L,
    nomad.opts=list(MAX_BB_EVAL=1L)),
    "conditional bandwidth search stopped: ambiguous numerical rank; row 14;")
  expect_identical(penalty(), before)
})
