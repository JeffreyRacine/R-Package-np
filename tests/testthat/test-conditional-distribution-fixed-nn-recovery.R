test_that("conditional CDF fixed MADS recovers automatic NN starts with the right grid caps", {
  expect_true(spawn_mpi_slaves(1L))
  on.exit(close_mpi_slaves(force=TRUE),add=TRUE)
  old <- options(np.messages=FALSE,np.tree=FALSE,np.extendednn=FALSE)
  on.exit(options(old),add=TRUE)
  x <- data.frame(x1=c(rep(0,32),sin((1:16)*sqrt(2))+(1:16)/32),
    x2=c(rep(0,32),sin((1:16)*sqrt(3))+(1:16)/48))
  y <- data.frame(y=sin((1:48)/3)+(1:48)/90)
  for(grid in c("generated","training","external"))for(type in c("generalized_nn","adaptive_nn")) {
    controls <- list(xdat=x,ydat=y,bwtype=type,bwmethod="cv.ls",regtype="lp",
      degree=c(1L,1L),nomad=FALSE,bwsolver="mads",nmulti=1L,itmax=20L,
      powell.remin=FALSE,ngrid=7L,nomad.opts=list(MAX_BB_EVAL=30L))
    if(grid=="training")controls$do.full.integral <- TRUE
    if(grid=="external")controls$gydat <- data.frame(y=seq(-1,2,length.out=7L)+.017)
    bw <- do.call(npcdistbw,controls)
    rr <- bw$nomad.restart.results
    expect_length(rr,2L)
    expect_true(isTRUE(rr[[2L]]$recovery))
    expect_identical(rr[[2L]]$recovery_witness$evaluations,3L)
    caps <- if(type=="adaptive_nn")rep(46,3L) else c(if(grid=="training")46 else 47,47,47)
    expect_equal(rr[[2L]]$start,caps,tolerance=0)
    expect_identical(bw$nomad.best.restart,2L)
    raw <- .npcdistbw_eval_only(xdat=x,ydat=y,gydat=controls$gydat,bws=bw,
      do.full.integral=isTRUE(controls$do.full.integral),ngrid=7L,invalid.penalty="dbmax")
    expect_true(.np_nn_raw_objective_valid(raw$objective))
    expect_equal(bw$fval,raw$objective,tolerance=2e-12)
    if(grid=="generated") {
      controls$bws <- rep(7,3L)
      controls$nomad.opts <- list(MAX_BB_EVAL=1L)
      expect_error(do.call(npcdistbw,controls),"did not return a raw-valid solution")
      controls$bws <- caps
      explicit <- do.call(npcdistbw,controls)
      expect_false(any(vapply(explicit$nomad.restart.results,function(z)isTRUE(z$recovery),logical(1L))))
      expect_true(.np_nn_raw_objective_valid(explicit$fval))
    }
  }
})

test_that("conditional CDF recovery caps retain response grid precedence and coordinate order", {
  template <- list(type="generalized_nn",iycon=c(TRUE,FALSE,TRUE))
  setup <- list(nobs=10L,cont_flat=c(1L,3L,5L))
  expect_identical(.npcdistbw_nn_recovery_caps(template,setup),rep(9L,3L))
  expect_identical(.npcdistbw_nn_recovery_caps(template,setup,do.full.integral=TRUE),c(8L,8L,9L))
  expect_identical(.npcdistbw_nn_recovery_caps(template,setup,gydat=data.frame(y=1),
    do.full.integral=TRUE),rep(9L,3L))
  template$type <- "adaptive_nn"
  expect_identical(.npcdistbw_nn_recovery_caps(template,setup,do.full.integral=TRUE),rep(8L,3L))
})
