# One tiny fixture covers the estimator-owned activity surface, not child totals.
npindex_activity_clock <- function() {
  state <- new.env(parent = emptyenv())
  state$time <- 0
  function() { state$time <- state$time + .25; state$time }
}

npindex_activity_fixture <- function() {
  set.seed(20260404)
  x <- data.frame(x1 = runif(24, -1, 1), x2 = runif(24, -1, 1))
  z <- x$x1 + .6 * x$x2
  list(x = x, y = sin(z) + .2*z^2 + rnorm(24, sd = .05))
}

test_that("index fitting owns initial activity, inference, prediction and bootstrap", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old <- options(np.messages = FALSE, np.tree = FALSE,
                 np.progress.interval.sec = 0)
  on.exit(options(old), add = TRUE)
  f <- npindex_activity_fixture()
  for (reg in c("lc", "lp")) {
    args <- list(xdat=f$x, ydat=f$y, bws=c(1,.6,.35),
                 method="ichimura", regtype=reg, bandwidth.compute=FALSE)
    if (reg=="lp") args$degree <- 1L
    bw <- do.call(npindexbw,args)
    for (uncertainty in c(FALSE,TRUE)) {
      call <- list(bws=bw,txdat=f$x,tydat=f$y,se=uncertainty)
      options(np.messages=FALSE)
      quiet <- do.call(npindex,call)
      options(np.messages=TRUE)
      actual <- capture_progress_shadow_trace(do.call(npindex,call),
        force_renderer="single_line",now=npindex_activity_clock())
      lines <- vapply(actual$trace,`[[`,"","line")
      events <- vapply(actual$trace,`[[`,"","event")
      expect_match(lines[1L],"Fitting single-index model.*elapsed 0.0s: fitted values")
      expect_false(any(grepl("Fitting regression",lines,fixed=TRUE)))
      expect_equal(sum(events=="finish"),1L)
      expect_equal(actual$value$mean,quiet$mean,tolerance=1e-14)
      expect_identical(actual$value$betavcov,quiet$betavcov)
      if(uncertainty)
        expect_true(any(grepl("coefficient covariance",lines,fixed=TRUE)))
    }
    pred <- capture_progress_shadow_trace(
      predict(quiet,newdata=f$x[c(2L,7L),],se.fit=TRUE),
      force_renderer="single_line",now=npindex_activity_clock())
    expect_match(pred$trace[[1L]]$line,"Fitting single-index model.*fitted values")
    expect_equal(sum(vapply(pred$trace,`[[`,"","event")=="finish"),1L)
  }
  actual <- capture_progress_shadow_trace(
    npindex(bws=bw,txdat=f$x,tydat=f$y,se.type="bootstrap",B=2L),
    force_renderer="single_line",now=npindex_activity_clock())
  lines <- vapply(actual$trace,`[[`,"","line")
  events <- vapply(actual$trace,`[[`,"","event")
  fit.end <- which(grepl("Fitting single-index model",lines,fixed=TRUE)&events=="finish")
  boot.start <- which(grepl("Bootstrapping single-index fit",lines,fixed=TRUE))[1L]
  expect_lt(fit.end,boot.start)
  expect_false(any(grepl("Fitting regression",lines,fixed=TRUE)))
  runtime <- getFromNamespace(".np_progress_runtime","npRmpi")
  expect_null(runtime$fit_forward)
  expect_null(runtime$fit_state)
})

test_that("search progress hands off once to the index activity owner", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  skip_if_not_installed("crs")
  old <- options(np.messages=TRUE,np.tree=FALSE,np.progress.interval.sec=0)
  on.exit(options(old),add=TRUE)
  f <- npindex_activity_fixture()
  dat <- cbind(f$x,y=f$y)
  for (nomad in c(FALSE,TRUE)) {
    actual <- capture_progress_shadow_trace(
      do.call(npindex,list(bws=y~x1+x2,data=dat,method="ichimura",
        nomad=nomad,degree.min=0L,degree.max=1L,nmulti=1L,se=FALSE)),
      force_renderer="single_line",now=npindex_activity_clock())
    lines <- vapply(actual$trace,`[[`,"","line")
    events <- vapply(actual$trace,`[[`,"","event")
    fit <- grep("Fitting single-index model",lines,fixed=TRUE)
    expect_gt(min(fit),1L)
    expect_true(any(grepl("bandwidth|Bandwidth|degree|Degree",lines[seq_len(min(fit)-1L)])))
    expect_equal(sum(events[fit]=="finish"),1L)
    expect_false(any(grepl("Fitting regression",lines[fit],fixed=TRUE)))
  }
})
