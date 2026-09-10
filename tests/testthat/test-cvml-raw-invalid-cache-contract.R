test_that("CVML raw invalid cache cannot become finite likelihood guidance", {
  ns <- asNamespace("np")
  old <- options(np.messages=FALSE, np.tree=FALSE, np.largeh=FALSE,
    np.objective.cache=TRUE, npRmpi.autodispatch.disable=TRUE)
  on.exit(options(old), add=TRUE)
  set.seed(1); x <- data.frame(x=runif(80L)); y <- data.frame(y=rbeta(80L,1,1))
  bw <- get("npcdensbw",ns)(xdat=x,ydat=y,regtype="lp",degree=3L,
    bws=c(.15,.12),bandwidth.compute=FALSE,cxkerbound="range",cykerbound="range")
  prep <- get(".npcdensbw_prepared_prepare_args",ns)(xdat=x,ydat=y,bws=bw,
    invalid.penalty="baseline",degree.search=TRUE)
  names(prep)[names(prep)=="penalty_mode"] <- "penalty.mode"
  names(prep)[names(prep)=="penalty_multiplier"] <- "penalty.multiplier"
  prepare <- get("npPreparedObjectivePrepareConditionalDensity",ns)
  guided <- get("npPreparedObjectiveEvalConditionalDensity",ns)
  raw <- get("npPreparedObjectiveEvalConditionalDensityRaw",ns)
  destroy <- get("npPreparedObjectiveDestroyConditionalDensity",ns)
  expect_true(do.call(prepare,prep))
  out <- tryCatch({
    bad <- c(.12,.15) # Native order X then Y, unlike public Y then X.
    g1 <- guided(bad,3L)
    r1 <- raw(bad,3L)
    r2 <- raw(bad,3L)
    healthy <- raw(c(.5,.5),0L)
    g2 <- guided(bad,3L)
    r3 <- raw(bad,3L)
    list(g1=g1,r1=r1,r2=r2,healthy=healthy,g2=g2,r3=r3)
  },finally=destroy())
  for(nm in c("r1","r2","r3")) {
    expect_identical(out[[nm]][[1L]],-.Machine$double.xmax)
    expect_length(out[[nm]],4L)
  }
  expect_identical(out$g1[[1L]],out$g2[[1L]])
  expect_true(is.finite(out$g1[[1L]]))
  expect_lt(abs(out$g1[[1L]]),.Machine$double.xmax)
  expect_true(is.finite(out$healthy[[1L]]))
  expect_lt(abs(out$healthy[[1L]]),.Machine$double.xmax)
  # Cached invalid evaluations must not execute the contribution guard again.
  expect_identical(out$r1[[4L]],0)
  expect_identical(out$r2[[4L]],0)
  expect_identical(out$r3[[4L]],0)
  expect_gt(out$g1[[4L]],0)
})
