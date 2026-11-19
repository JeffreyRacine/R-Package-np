test_that("pairs plot layouts acquire locally while wild retains its owner", {
  if (!spawn_mpi_slaves(n=1L)) skip("MPI pool unavailable in this test context")
  on.exit(close_mpi_slaves(force=TRUE), add=TRUE)
  old <- options(np.messages=FALSE, np.plot.progress=FALSE,
    npRmpi.local.regression.mode=FALSE, npRmpi.autodispatch.disable=FALSE)
  on.exit(options(old), add=TRUE)
  i <- seq_len(40L)
  x <- data.frame(x=sin(i*sqrt(2)),z=cos(i*sqrt(3)))
  y <- x$x+.4*x$z^2+.1*sin(i*sqrt(5))
  bw <- npregbw(xdat=x,ydat=y,bws=c(.4,.6),bandwidth.compute=FALSE,
    regtype="lp",degree=c(1L,2L),bernstein.basis=FALSE)
  fit <- npreg(bws=bw,txdat=x,tydat=y)
  real.direct <- npRmpi:::.np_regression_direct
  captured <- public <- 0L
  real.public <- npRmpi:::npreg
  testthat::local_mocked_bindings(
    npreg=function(...) { public <<- public+1L; real.public(...) },
    .np_regression_direct=function(...) {
      args <- list(...)
      if (isTRUE(args$fit.inference)) {
        captured <<- captured+1L
        expect_true(isTRUE(getOption("npRmpi.local.regression.mode")))
        expect_true(isTRUE(getOption("npRmpi.autodispatch.disable")))
      }
      real.direct(...)
    }, .package="npRmpi")
  for (law in c("inid","wild")) for (layout in c("mean","gradient","surface")) {
    captured <- public <- 0L
    prior.options <- options("npRmpi.local.regression.mode", "npRmpi.autodispatch.disable")
    out <- plot(fit, errors="bootstrap", center="bias-corrected",
      bootstrap=law, B=39L, alpha=.2, neval=3L, band="pointwise",
      gradients=layout=="gradient", perspective=layout=="surface",
      view="fixed", output="data", random.seed=107L,
      bias.control=np_bias_control(adaptation="none"))
    expect_s3_class(out,"np_plot_data")
    expected <- if(layout=="surface")1L else 2L
    expect_identical(captured,if(law=="inid")expected else 0L)
    expect_identical(public,if(law=="wild")expected else 0L)
    expect_identical(do.call(options,as.list(names(prior.options))),prior.options)
  }
})
