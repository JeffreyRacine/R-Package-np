test_that("external CDF recovery separates deleted and incumbent count domains", {
  template <- list(type="generalized_nn",iycon=TRUE,
    cxkertype="gaussian",cykertype="gaussian")
  setup <- list(nobs=24L,cont_flat=1:2)
  expect_identical(.npcdistbw_nn_recovery_caps(template,setup),c(22L,22L))
  expect_identical(.npcdistbw_nn_recovery_caps(template,setup,incumbent=TRUE),c(23L,23L))
  for(k in 21:24) {
    expected <- if(k %in% c(21,23))list(c(22,22,.3))else list()
    expect_identical(.np_nn_ordinary_schedule(c(k,k,.3),1:2,c(22,22),c(23,23)),expected)
  }
  expect_identical(.npcdistbw_nn_recovery_caps(template,setup,
    do.full.integral=TRUE),c(22L,23L))
  for(role in c("cxkertype","cykertype")) {
    beta <- template; beta[[role]] <- "beta"
    expect_identical(.npcdistbw_nn_recovery_caps(beta,setup),c(23L,23L))
  }
  template$type <- "adaptive_nn"
  expect_identical(.npcdistbw_nn_recovery_caps(template,setup,incumbent=TRUE),c(22L,22L))
})

test_that("the external recovery endpoint is certified by literal deleted fits", {
  old <- options(np.messages=FALSE,np.tree=FALSE,np.extendednn=FALSE)
  on.exit(options(old),add=TRUE)
  x <- data.frame(x=c(rep(0,22),1,2))
  y <- data.frame(y=seq(.1,.9,length.out=24))
  grid <- data.frame(y=c(.15,.5,.85))
  b <- npcdistbw(xdat=x,ydat=y,bws=c(22,22),bwtype="generalized_nn",
    regtype="lc",bandwidth.compute=FALSE)
  raw <- .npcdistbw_eval_only(x,y,gydat=grid,bws=b,invalid.penalty="dbmax")$objective
  literal <- 0
  for(i in 1:24) {
    fit <- fitted(npcdist(bws=b,txdat=x[-i,,drop=FALSE],tydat=y[-i,,drop=FALSE],
      exdat=x[rep(i,3),,drop=FALSE],eydat=grid))
    literal <- literal+mean((as.numeric(y$y[i]<=grid$y)-fit)^2)/24
  }
  expect_equal(raw,literal,tolerance=2e-12)
  for(solver in c("powell","mads","mads+powell")) {
    set.seed(42)
    fit <- npcdistbw(xdat=x,ydat=y,gydat=grid,bwtype="generalized_nn",
      regtype="lc",bwsolver=solver,nmulti=1L,itmax=20L,powell.remin=FALSE,
      nomad.opts=list(MAX_BB_EVAL=30L))
    got <- .npcdistbw_eval_only(x,y,gydat=grid,bws=fit,invalid.penalty="dbmax")$objective
    expect_true(.np_nn_raw_objective_valid(got))
    expect_equal(fit$fval,got,tolerance=2e-12)
    if(solver!="powell")expect_true(any(vapply(fit$nomad.restart.results,
      function(z)isTRUE(z$recovery),logical(1))))
  }
})
