test_that("surface overlays respect explicit height limits and retain statistics", {
  if (!spawn_mpi_slaves(n=1L)) skip("MPI pool unavailable in this test context")
  on.exit(close_mpi_slaves(force=TRUE),add=TRUE)
  old <- options(np.messages=FALSE,plot.par.mfrow=FALSE)
  on.exit(options(old),add=TRUE)
  grDevices::pdf(NULL); on.exit(grDevices::dev.off(),add=TRUE)
  received <- list()
  range.calls <- 0L
  original.range <- .np_plot_overlay_range
  testthat::local_mocked_bindings(.np_plot_overlay_range=function(...) {
    range.calls <<- range.calls + 1L
    original.range(...)
  }, .package="npRmpi")
  original.persp <- getS3method("persp","default",envir=asNamespace("graphics"))
  testthat::local_mocked_bindings(persp.default=function(x,y,z,zlim,...) {
    received[[length(received)+1L]] <<- zlim
    original.persp(x=x,y=y,z=z,zlim=zlim,...)
  },.package="graphics")
  d <- data.frame(x=seq(.05,.95,length.out=30L),z=cos(seq(.1,5,length.out=30L)))
  d$y <- sin(3*d$x)+.2*d$z
  fits <- list(npreg(y~x+z,data=d,bws=c(.3,.4),bandwidth.compute=FALSE),
    npplreg(y~x|z,data=d,bws=matrix(c(.4,.4),2L,1L),bandwidth.compute=FALSE),
    npscoef(y~x|z,data=d,bws=.4,bandwidth.compute=FALSE))
  payload <- function(x) lapply(x$panels,function(p)p[c("grid","variables","values","intervals")])
  for(fit in fits) {
    original <- fit
    received <- list()
    range.calls <- 0L
    plot(fit, neval=4L, perspective=TRUE, view="fixed", data.overlay=TRUE)
    expect_true(length(received) >= 2L)
    expect_gt(range.calls, 0L)
    expect_true(all(vapply(received, function(limits)
      limits[1L] <= min(d$y) && limits[2L] >= max(d$y), logical(1L))))
    ordinary <- plot(fit,neval=4L,perspective=TRUE,view="fixed",data.overlay=TRUE,output="data")
    received <- list()
    expect_warning(limited <- plot(fit,neval=4L,perspective=TRUE,view="fixed",data.overlay=TRUE,
      zlim=c(.45,.6),output="plot-data"), "surface extends beyond the box")
    expect_true(length(received)>=2L)
    expect_true(all(vapply(received,identical,logical(1L),c(.45,.6))))
    expect_identical(payload(limited),payload(ordinary))
    expect_identical(fit,original)
    received <- list()
    expect_warning(plot(fit,neval=4L,perspective=TRUE,view="fixed",data.overlay=FALSE,zlim=c(.45,.6)),
      "surface extends beyond the box")
    expect_true(all(vapply(received,identical,logical(1L),c(.45,.6))))
  }
})
