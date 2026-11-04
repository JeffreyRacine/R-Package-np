test_that("continuous regression panels forward xlim without changing inference", {
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
                 bias_control=np_bias_control(degree=2),gaussian_draws=31L,
                 neval=3L,common_scale=common,gradients=grad)
    reference <- do.call(plot,c(args,list(output="data")))
    seed <- .Random.seed
    value <- do.call(plot,c(args,list(output="plot-data",xlim=c(.8,.2))))
    expect_identical(received,c(.8,.2))
    expect_identical(value$inference, reference$inference)
    expect_identical(value$targets, reference$targets)
    for (i in seq_along(value$panels))
      expect_identical(value$panels[[i]][c("grid","variables","values","intervals")],
                       reference$panels[[i]][c("grid","variables","values","intervals")])
    expect_identical(value$panels[[1L]]$display$style$xlim, c(.8,.2))
    expect_identical(.Random.seed,seed)
  }
})
