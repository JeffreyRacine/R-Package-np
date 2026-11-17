test_that("regression capture compares model state independently of call provenance", {
  same <- npRmpi:::.np_plot_regression_capture_same_bws
  a <- list(bw=c(.5,.3),degree=c(1L,2L),type="fixed",call=quote(npregbw(y~x)))
  attr(a$call,".Environment") <- environment()
  b <- a;attr(b$call,".Environment") <- NULL
  expect_true(same(a,b))
  before <- a;expect_true(same(a,a));expect_identical(a,before)
  for(field in c("bw","degree","type")) {
    changed <- b;changed[[field]] <- switch(field,bw=c(.5,.4),degree=c(1L,3L),type="adaptive_nn")
    expect_false(same(a,changed))
  }
  tagged <- b;attr(tagged,"npRmpi.dispatch.mode") <- "auto"
  attr(tagged$bw,"npRmpi.autodispatch.remote") <- "transport-only"
  expect_true(same(a,tagged))
  tagged <- b;attr(tagged,"unknown.model.state") <- 1L
  expect_false(same(a,tagged))
  b$new.computational.field <- TRUE
  expect_false(same(a,b))
  expect_false(same(NULL,a))
})


test_that("capture identity excludes only top-level diagnostic profile", {
  same <- npRmpi:::.np_plot_regression_capture_same_bws
  a <- list(bw=.5, degree=1L, type="fixed", nested=list(timing.profile=1))
  before <- a
  b <- a; b$timing.profile <- list(total=10, where="npregbw")
  expect_true(same(a,b))
  expect_true(same(b,a))
  expect_identical(a,before)
  saved <- b
  expect_true(same(b,b))
  expect_identical(b,saved)
  b$timing.profile <- list(total=20, where="npreg")
  expect_true(same(saved,b))
  for (field in names(a)) {
    bad <- b; bad[[field]] <- list(changed=bad[[field]])
    expect_false(same(a,bad), info=field)
  }
  bad <- b; bad$timing <- 1
  expect_false(same(a,bad))
  bad <- b; bad$nested$timing.profile <- 2
  expect_false(same(a,bad))
  bad <- b; attr(bad,"timing.profile") <- 1
  expect_false(same(a,bad))
  bad <- b; bad$new.model.field <- 1
  expect_false(same(a,bad))
  bad <- b; class(bad) <- "different.bandwidth.class"
  expect_false(same(a,bad))
})

test_that("native bandwidth plots preserve captured inference with telemetry", {
  if (!spawn_mpi_slaves(n=1L)) skip("MPI pool unavailable in this test context")
  on.exit(close_mpi_slaves(force=TRUE),add=TRUE)
  old <- options(np.messages=FALSE, np.tree=FALSE)
  on.exit(options(old),add=TRUE)
  set.seed(42)
  x <- data.frame(x=runif(60,-2,2)); y <- x$x+rnorm(60)
  for (scale in c(FALSE,TRUE)) {
    bw <- npregbw(xdat=x, ydat=y, bws=if(scale)10000 else .5,
      bwscaling=scale, regtype="ll", bandwidth.compute=FALSE)
    saved <- bw
    fitted.reference <- npreg(bw)
    seed <- .Random.seed
    for (center in c("estimate","bias-corrected")) {
      value <- plot(bw, errors="asymptotic", center=center, band="all",
        neval=3L, gaussian.draws=63L, output="data")
      expect_s3_class(value,"np_plot_data")
      reference <- plot(fitted.reference, errors="asymptotic", center=center,
        band="all", neval=3L, gaussian.draws=63L, output="data")
      for (selected.band in c("pointwise","simultaneous","bonferroni"))
        expect_identical(as.data.frame(value, panel=1L, band=selected.band),
          as.data.frame(reference, panel=1L, band=selected.band))
      expect_identical(.Random.seed,seed)
      expect_identical(bw,saved)
    }
  }
})
