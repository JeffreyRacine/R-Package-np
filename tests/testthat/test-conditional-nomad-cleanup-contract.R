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
  # Pin the original failing point: local suite mode can take a different
  # optimizer trajectory. Distributed free-degree recovery is proved separately.
  for (visible in c(FALSE, TRUE)) {
    options(np.messages=visible)
    set.seed(20260322)
    dat <- data.frame(x=sort(runif(14)), y=sort(runif(14)))
    expect_error(npcdensbw(y~x, data=dat, regtype='lp', degree=1L,
      bws=c(0.048130207494902745,0.025634352719319447),
      bwsolver='mads', nmulti=1L, nomad.opts=list(MAX_BB_EVAL=8L)),
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


test_that("conditional CDF block failure shares numerical metadata", {
  source_file <- file.path(npRmpi_namespace_hygiene_root(), "src", "jksum.c")
  source <- paste(readLines(source_file, warn=FALSE), collapse="\n")
  start <- regexpr("np_conditional_distribution_cvls_lp_one\\(", source)[1L]
  end <- regexpr("#undef NP_CDIST_ONEBLOCK_ALIGN", source, fixed=TRUE)[1L]
  expect_gt(start, 0L)
  expect_gt(end, start)
  owner <- substr(source, start, end-1L)
  expect_match(owner, "np_conditional_failure_reduce(1, local_fail)", fixed=TRUE)
  expect_false(grepl("MPI_Allreduce(&local_fail", owner, fixed=TRUE))
})
