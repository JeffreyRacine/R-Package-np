test_that("continuous regression panels forward xlim without changing inference", {
  if (!spawn_mpi_slaves(n=1L)) skip("MPI pool unavailable in this test context")
  on.exit(close_mpi_slaves(force=TRUE), add=TRUE)
  old <- options(np.messages=FALSE, plot.par.mfrow=FALSE)
  on.exit(options(old), add=TRUE)
  grDevices::pdf(tempfile(fileext=".pdf"))
  on.exit(grDevices::dev.off(), add=TRUE)
  received <- NULL
  original.window <- graphics::plot.window
  testthat::local_mocked_bindings(
    plot.window=function(xlim, ...) {
      received <<- xlim
      original.window(xlim=xlim, ...)
    }, .package="graphics")
  set.seed(49620)
  x <- data.frame(x=runif(40)); y <- sin(x$x*3)+rnorm(40,sd=.2)
  bw <- npregbw(xdat=x,ydat=y,bws=.3,regtype="ll",bandwidth.compute=FALSE)
  fit <- npreg(bw,txdat=x,tydat=y)
  for (common in c(FALSE,TRUE)) for (grad in c(FALSE,TRUE)) {
    args <- list(x=fit,errors="asymptotic",band="all",center="bias-corrected",
                 # Keep this renderer test on the qualified multi-chunk route.
                 # The inherited single-chunk MPI argument collision is tracked
                 # separately; changing its dispatcher is outside this repair.
                 bias.control=np_bias_control(degree=2),gaussian.draws=63L,
                 neval=3L,common.scale=common,gradients=grad)
    reference <- do.call(plot,c(args,list(output="data")))
    seed <- .Random.seed
    value <- do.call(plot,c(args,list(output="plot-data",xlim=c(.8,.2))))
    expect_identical(received,c(.8,.2))
    expect_identical(value$targets, reference$targets)
    expect_identical(value$inference, reference$inference)
    for (j in seq_along(value$panels)) {
      expect_identical(value$panels[[j]]$grid, reference$panels[[j]]$grid)
      expect_identical(value$panels[[j]]$values, reference$panels[[j]]$values)
      expect_identical(value$panels[[j]]$intervals, reference$panels[[j]]$intervals)
    }
    expect_identical(.Random.seed,seed)
  }
})
