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
                 bias.control=np_bias_control(degree=2),gaussian.draws=31L,
                 neval=3L,common.scale=common,gradients=grad)
    expect_warning(
      reference <- do.call(plot,c(args,list(output="data"))),
      "Empirical band quantiles have limited draw resolution (bonferroni)",
      fixed = TRUE)
    seed <- .Random.seed
    expect_warning(
      value <- do.call(plot,c(args,list(output="plot-data",xlim=c(.8,.2)))),
      "Empirical band quantiles have limited draw resolution (bonferroni)",
      fixed = TRUE)
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

test_that("regression surface wires retain coordinates and asymmetric widths", {
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  group <- ordered(rep(1:2, each = 24))
  x <- rep(seq(0, 1, length.out = 24), 2)
  fit <- npreg(y ~ group + x, data = data.frame(group, x, y = as.integer(group) + x),
               bws = c(.2, .3), bandwidth.compute = FALSE)
  seen <- NULL
  widths <- cbind(1:6, 6:1)
  scopes <- list(pointwise = widths, simultaneous = 2 * widths, bonferroni = 3 * widths)
  testthat::local_mocked_bindings(
    .np_plot_regression_eval = function(...) list(mean = 10 * (1:6), merr = numeric(6)),
    compute.bootstrap.errors = function(...) list(
      boot.err = cbind(widths, NA_real_), boot.all.err = scopes),
    .np_plot_render_surface_rgl = function(draw.extras, ...) {
      draw.extras()
      invisible(NULL)
    },
    .np_plot_error_surfaces_rgl = function(x, y, lerr, herr, lerr.all, herr.all, ...) {
      seen <<- list(x = x, y = y, lower = lerr, upper = herr,
                    lower.all = lerr.all, upper.all = herr.all)
    }, .package = "np")
  plot(fit, renderer = "rgl", errors = "bootstrap", band = "all",
       neval = 3L, B = 19L, data.rug = FALSE, data.overlay = FALSE)
  expect_length(seen$x, 2L)
  expect_length(seen$y, 3L)
  expect_identical(seen$lower, matrix(c(9, 18, 27, 36, 45, 54), 2, 3))
  expect_identical(seen$upper, matrix(c(16, 25, 34, 43, 52, 61), 2, 3))
  expect_named(seen$lower.all, names(scopes))
  expect_named(seen$upper.all, names(scopes))
  for (name in names(scopes)) {
    expect_identical(seen$lower.all[[name]], matrix(10 * (1:6) - scopes[[name]][, 1], 2))
    expect_identical(seen$upper.all[[name]], matrix(10 * (1:6) + scopes[[name]][, 2], 2))
  }
})
