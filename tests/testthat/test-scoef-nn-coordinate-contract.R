test_that("NN smoothing maps preserve categorical coordinates and their bounds", {
  expect_true(spawn_mpi_slaves(1L))
  on.exit(close_mpi_slaves(force=TRUE),add=TRUE)
  old <- options(np.extendednn=FALSE,np.messages=FALSE)
  on.exit(options(old),add=TRUE)
  for(type in c("generalized_nn","adaptive_nn")) {
    for(order in list(1:3,c(2L,1L,3L),c(1L,3L,2L))) {
      icon <- c(FALSE,TRUE,FALSE)[order]
      param <- c(.2,11.6,.3)[order]
      expected <- c(.2,12,.3)[order]
      upper <- ifelse(icon,Inf,.5)
      expect_identical(.npscoef_nn_candidate_bandwidth(param,type,80L,icon),expected)
      expect_identical(.np_scbandwidth_manual_nn_validate(expected,80L,icon),expected)
      sc <- list(type=type,bw=expected,icon=icon,nobs=80L,scaling=FALSE,bandwidth=list(expected))
      applied <- .npscoef_apply_bw_to_scbw(sc,param)
      expect_identical(applied$bw,expected)
      expect_identical(applied$bandwidth[[1L]],expected)
      expect_true(.npscoef_candidate_is_admissible(expected,type,80L,upper=upper,icon=icon))
      bad <- expected;bad[which(!icon)[1L]] <- .6
      expect_false(.npscoef_candidate_is_admissible(bad,type,80L,upper=upper,icon=icon))
      bad[which(!icon)[1L]] <- -.1
      expect_false(.npscoef_candidate_is_admissible(bad,type,80L,upper=upper,icon=icon))
    }
    expect_error(.npscoef_nn_candidate_bandwidth(c(.2,12),type,80L,TRUE),"invalid nearest-neighbor coordinate map")
    expect_error(.np_scbandwidth_manual_nn_validate(c(.2,12),80L,TRUE),"invalid nearest-neighbor coordinate map")
    expect_identical(.npscoef_nn_candidate_bandwidth(c(1,1.5,2.5,3.5,80),type,80L),
      c(NA_real_,2,2,4,NA_real_))
    expect_identical(.npscoef_nn_candidate_bandwidth(c(.2,.3),type,80L,c(FALSE,FALSE)),c(.2,.3))
  }
})

test_that("NN start helpers preserve old continuous RNG and typed categorical draws", {
  old <- options(np.extendednn=FALSE)
  on.exit(options(old),add=TRUE)
  for(type in c("fixed","generalized_nn","adaptive_nn")) {
    args <- list(param=rep(8,3L),bwtype=type,nobs=80L)
    mapped <- c(args,list(icon=rep(TRUE,3L),iord=rep(FALSE,3L),iuno=rep(FALSE,3L)))
    expect_identical(do.call(.npscoef_default_start_bandwidth,args),do.call(.npscoef_default_start_bandwidth,mapped))
    set.seed(9303L);a <- do.call(.npscoef_random_start_bandwidth,args);after <- .Random.seed
    set.seed(9303L);b <- do.call(.npscoef_random_start_bandwidth,mapped)
    expect_identical(a,b)
    expect_identical(after,.Random.seed)
  }
  for(type in c("generalized_nn","adaptive_nn")) {
    args <- list(param=c(.25,8,.4),bwtype=type,nobs=80L,icon=c(FALSE,TRUE,FALSE),
      iuno=c(TRUE,FALSE,FALSE),iord=c(FALSE,FALSE,TRUE))
    expect_identical(do.call(.npscoef_default_start_bandwidth,args),c(.25,9,.4))
    set.seed(9303L);a <- do.call(.npscoef_random_start_bandwidth,args);after <- .Random.seed
    set.seed(9303L);k <- round(runif(1L,2,79));cats <- c(.25,.4)*runif(2L,.5,1.5)
    expect_identical(a,c(cats[1L],k,cats[2L]))
    expect_identical(after,.Random.seed)
  }
})

test_that("semihat preserves valid full-width bounds and the scalar index case", {
  z <- data.frame(u=factor(c("a","b","a","b")),x=c(-1,0,1,2),
    o=ordered(c("lo","hi","lo","hi")))
  source <- list(regtype="lp",basis="glp",degree=1L,bernstein.basis=FALSE,
    type="fixed",ckerbound="fixed",ckerlb=c(-Inf,-2,-Inf),ckerub=c(Inf,3,Inf))
  args <- .np_semihat_make_regbw_args(source,z,1:4,c(.2,.5,.3))
  expect_identical(args$ckerlb,source$ckerlb)
  expect_identical(args$ckerub,source$ckerub)
  source$ckerlb <- -2;source$ckerub <- 3
  args <- .np_semihat_make_regbw_args(source,z["x"],1:4,.5)
  expect_identical(args$ckerlb,-2)
  expect_identical(args$ckerub,3)
})
