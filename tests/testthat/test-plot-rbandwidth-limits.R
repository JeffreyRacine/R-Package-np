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

test_that("coordinator draws the published regression surface bounds", {
  if (!spawn_mpi_slaves(n = 1L)) skip("MPI pool unavailable in this test context")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  xy <- expand.grid(x = seq(-1, 1, length.out = 7L), z = seq(0, 1, length.out = 8L))
  xy$y <- with(xy, x^2 + z + sin(5 * x + z))
  fit <- npreg(y ~ x + z, data = xy, bws = c(.3, .4), bandwidth.compute = FALSE)
  captured <- NULL
  testthat::local_mocked_bindings(
    .np_plot_render_surface_rgl = function(draw.extras, ...) {
      draw.extras()
      invisible(NULL)
    },
    .np_plot_error_surfaces_rgl = function(lerr, herr, ...) {
      captured <<- list(lower = lerr, upper = herr)
    }, .package = "npRmpi")
  out <- plot(fit, renderer = "rgl", output = "plot-data", errors = "asymptotic",
              band = "pointwise", neval = 4L, data.rug = FALSE, data.overlay = FALSE)
  intervals <- out$panels[[1L]]$intervals
  expect_identical(dim(captured$lower), c(4L, 4L))
  expect_identical(as.vector(captured$lower), intervals$lower)
  expect_identical(as.vector(captured$upper), intervals$upper)
})
