quiet_capture <- function(expr) {
  out <- NULL
  sinkfile <- tempfile()
  on.exit(unlink(sinkfile), add = TRUE)
  invisible(capture.output(out <- eval.parent(substitute(expr)), file = sinkfile))
  out
}

flatten_plot_mean <- function(plot_obj) {
  vals <- unlist(lapply(plot_obj, function(comp) {
    keys <- c("mean", "dens", "dist", "condens", "condist")
    hit <- keys[keys %in% names(comp)]
    if (!length(hit))
      return(numeric(0))
    as.double(comp[[hit[1L]]])
  }), use.names = FALSE)
  vals[is.finite(vals)]
}

run_bootstrap_plot <- function(bw,
                               xdat = NULL,
                               ydat = NULL,
                               zdat = NULL,
                               boot_num = 9L,
                               boot_method = "wild") {
  args <- list(
    bw,
    plot.behavior = "data",
    perspective = FALSE,
    plot.errors.method = "bootstrap",
    plot.errors.boot.method = boot_method,
    plot.errors.boot.num = boot_num
  )
  if (!is.null(xdat))
    args$xdat <- xdat
  if (!is.null(ydat))
    args$ydat <- ydat
  if (!is.null(zdat))
    args$zdat <- zdat

  quiet_capture(suppressWarnings(do.call(plot, args)))
}

test_that("rbandwidth plot helper uses method-defining kernel options", {
  set.seed(9301)
  n <- 50L
  x <- runif(n)
  y <- sin(2 * pi * x) + rnorm(n, sd = 0.1)
  xdat <- data.frame(x = x)

  bw.gaussian <- npregbw(
    y ~ x,
    data = data.frame(y = y, x = x),
    bws = 0.2,
    bandwidth.compute = FALSE,
    regtype = "lp",
    basis = "glp",
    degree = 2L,
    bernstein.basis = TRUE,
    bwtype = "fixed",
    ckertype = "gaussian"
  )

  bw.epan <- npregbw(
    y ~ x,
    data = data.frame(y = y, x = x),
    bws = 0.2,
    bandwidth.compute = FALSE,
    regtype = "lp",
    basis = "glp",
    degree = 2L,
    bernstein.basis = TRUE,
    bwtype = "fixed",
    ckertype = "epanechnikov",
    ckerorder = 2L
  )

  set.seed(9302)
  out.gaussian <- run_bootstrap_plot(bw.gaussian, xdat = xdat, ydat = y, boot_num = 9L)
  set.seed(9302)
  out.epan <- run_bootstrap_plot(bw.epan, xdat = xdat, ydat = y, boot_num = 9L)

  mean.gaussian <- flatten_plot_mean(out.gaussian)
  mean.epan <- flatten_plot_mean(out.epan)

  expect_gt(length(mean.gaussian), 0L)
  expect_equal(length(mean.gaussian), length(mean.epan))
  expect_gt(max(abs(mean.gaussian - mean.epan)), 1e-6)
})

test_that("sibandwidth plot helper uses method-defining kernel options", {
  set.seed(9303)
  n <- 60L
  xdat <- data.frame(x1 = runif(n), x2 = runif(n))
  y <- xdat$x1^2 + 0.5 * xdat$x2 + rnorm(n, sd = 0.1)

  bw.gaussian <- npindexbw(
    xdat = xdat,
    ydat = y,
    bws = c(0.2, 0.2, 0.25),
    bandwidth.compute = FALSE,
    regtype = "lp",
    basis = "glp",
    degree = 2L,
    bernstein.basis = TRUE,
    bwtype = "fixed",
    ckertype = "gaussian"
  )

  bw.epan <- npindexbw(
    xdat = xdat,
    ydat = y,
    bws = c(0.2, 0.2, 0.25),
    bandwidth.compute = FALSE,
    regtype = "lp",
    basis = "glp",
    degree = 2L,
    bernstein.basis = TRUE,
    bwtype = "fixed",
    ckertype = "epanechnikov",
    ckerorder = 2L
  )

  set.seed(9304)
  out.gaussian <- run_bootstrap_plot(bw.gaussian, xdat = xdat, ydat = y, boot_num = 9L)
  set.seed(9304)
  out.epan <- run_bootstrap_plot(bw.epan, xdat = xdat, ydat = y, boot_num = 9L)

  mean.gaussian <- flatten_plot_mean(out.gaussian)
  mean.epan <- flatten_plot_mean(out.epan)

  expect_gt(length(mean.gaussian), 0L)
  expect_equal(length(mean.gaussian), length(mean.epan))
  expect_gt(max(abs(mean.gaussian - mean.epan)), 1e-6)
})

test_that("scbandwidth plot helper uses method-defining kernel options", {
  set.seed(9305)
  n <- 55L
  xdat <- data.frame(x = runif(n))
  zdat <- data.frame(z = runif(n))
  y <- sin(2 * pi * xdat$x) * zdat$z + rnorm(n, sd = 0.1)

  bw.gaussian <- npscoefbw(
    xdat = xdat,
    zdat = zdat,
    ydat = y,
    bws = 0.2,
    bandwidth.compute = FALSE,
    regtype = "lp",
    basis = "glp",
    degree = 2L,
    bernstein.basis = TRUE,
    bwtype = "fixed",
    ckertype = "gaussian"
  )

  bw.epan <- npscoefbw(
    xdat = xdat,
    zdat = zdat,
    ydat = y,
    bws = 0.2,
    bandwidth.compute = FALSE,
    regtype = "lp",
    basis = "glp",
    degree = 2L,
    bernstein.basis = TRUE,
    bwtype = "fixed",
    ckertype = "epanechnikov",
    ckerorder = 2L
  )

  set.seed(9306)
  out.gaussian <- run_bootstrap_plot(bw.gaussian, xdat = xdat, ydat = y, zdat = zdat, boot_num = 9L)
  set.seed(9306)
  out.epan <- run_bootstrap_plot(bw.epan, xdat = xdat, ydat = y, zdat = zdat, boot_num = 9L)

  mean.gaussian <- flatten_plot_mean(out.gaussian)
  mean.epan <- flatten_plot_mean(out.epan)

  expect_gt(length(mean.gaussian), 0L)
  expect_equal(length(mean.gaussian), length(mean.epan))
  expect_gt(max(abs(mean.gaussian - mean.epan)), 1e-6)
})

test_that("plbandwidth plot helper uses method-defining kernel options", {
  set.seed(9307)
  n <- 40L
  xdat <- data.frame(x = runif(n))
  zdat <- data.frame(z = runif(n))
  y <- 2 * xdat$x + sin(2 * pi * zdat$z) + rnorm(n, sd = 0.1)

  bw.gaussian <- quiet_capture(npplregbw(
    xdat = xdat,
    zdat = zdat,
    ydat = y,
    nmulti = 1L,
    regtype = "lp",
    basis = "glp",
    degree = 2L,
    bernstein.basis = TRUE,
    bwtype = "fixed",
    ckertype = "gaussian"
  ))

  bw.epan <- quiet_capture(npplregbw(
    xdat = xdat,
    zdat = zdat,
    ydat = y,
    nmulti = 1L,
    regtype = "lp",
    basis = "glp",
    degree = 2L,
    bernstein.basis = TRUE,
    bwtype = "fixed",
    ckertype = "epanechnikov",
    ckerorder = 2L
  ))

  set.seed(9308)
  out.gaussian <- run_bootstrap_plot(bw.gaussian, xdat = xdat, ydat = y, zdat = zdat, boot_num = 7L)
  set.seed(9308)
  out.epan <- run_bootstrap_plot(bw.epan, xdat = xdat, ydat = y, zdat = zdat, boot_num = 7L)

  mean.gaussian <- flatten_plot_mean(out.gaussian)
  mean.epan <- flatten_plot_mean(out.epan)

  expect_gt(length(mean.gaussian), 0L)
  expect_equal(length(mean.gaussian), length(mean.epan))
  expect_gt(max(abs(mean.gaussian - mean.epan)), 1e-6)
})

test_that("bandwidth plot helper uses method-defining kernel options", {
  set.seed(9309)
  n <- 70L
  xdat <- data.frame(x = runif(n))

  bw.gaussian <- npudensbw(
    dat = xdat,
    bws = 0.25,
    bandwidth.compute = FALSE,
    bwtype = "fixed",
    ckertype = "gaussian"
  )

  bw.epan <- npudensbw(
    dat = xdat,
    bws = 0.25,
    bandwidth.compute = FALSE,
    bwtype = "fixed",
    ckertype = "epanechnikov",
    ckerorder = 2L
  )

  set.seed(9310)
  out.gaussian <- run_bootstrap_plot(
    bw.gaussian,
    xdat = xdat,
    boot_num = 9L,
    boot_method = "inid"
  )
  set.seed(9310)
  out.epan <- run_bootstrap_plot(
    bw.epan,
    xdat = xdat,
    boot_num = 9L,
    boot_method = "inid"
  )

  mean.gaussian <- flatten_plot_mean(out.gaussian)
  mean.epan <- flatten_plot_mean(out.epan)

  expect_gt(length(mean.gaussian), 0L)
  expect_equal(length(mean.gaussian), length(mean.epan))
  expect_gt(max(abs(mean.gaussian - mean.epan)), 1e-6)
})

test_that("dbandwidth plot helper uses method-defining kernel options", {
  set.seed(9311)
  n <- 70L
  xdat <- data.frame(x = runif(n))

  bw.gaussian <- npudistbw(
    dat = xdat,
    bws = 0.25,
    bandwidth.compute = FALSE,
    bwtype = "fixed",
    ckertype = "gaussian"
  )

  bw.epan <- npudistbw(
    dat = xdat,
    bws = 0.25,
    bandwidth.compute = FALSE,
    bwtype = "fixed",
    ckertype = "epanechnikov",
    ckerorder = 2L
  )

  set.seed(9312)
  out.gaussian <- run_bootstrap_plot(
    bw.gaussian,
    xdat = xdat,
    boot_num = 9L,
    boot_method = "inid"
  )
  set.seed(9312)
  out.epan <- run_bootstrap_plot(
    bw.epan,
    xdat = xdat,
    boot_num = 9L,
    boot_method = "inid"
  )

  mean.gaussian <- flatten_plot_mean(out.gaussian)
  mean.epan <- flatten_plot_mean(out.epan)

  expect_gt(length(mean.gaussian), 0L)
  expect_equal(length(mean.gaussian), length(mean.epan))
  expect_gt(max(abs(mean.gaussian - mean.epan)), 1e-6)
})

test_that("conbandwidth plot helper uses method-defining x/y kernel options", {
  set.seed(9313)
  n <- 64L
  xdat <- data.frame(x = runif(n))
  ydat <- data.frame(y = runif(n))

  bw.gaussian <- npcdensbw(
    xdat = xdat,
    ydat = ydat,
    bws = c(0.22, 0.22),
    bandwidth.compute = FALSE,
    bwtype = "fixed",
    cxkertype = "gaussian",
    cykertype = "gaussian"
  )

  bw.epan <- npcdensbw(
    xdat = xdat,
    ydat = ydat,
    bws = c(0.22, 0.22),
    bandwidth.compute = FALSE,
    bwtype = "fixed",
    cxkertype = "epanechnikov",
    cykertype = "epanechnikov",
    cxkerorder = 2L,
    cykerorder = 2L
  )

  set.seed(9314)
  out.gaussian <- run_bootstrap_plot(
    bw.gaussian,
    xdat = xdat,
    ydat = ydat,
    boot_num = 9L,
    boot_method = "inid"
  )
  set.seed(9314)
  out.epan <- run_bootstrap_plot(
    bw.epan,
    xdat = xdat,
    ydat = ydat,
    boot_num = 9L,
    boot_method = "inid"
  )

  mean.gaussian <- flatten_plot_mean(out.gaussian)
  mean.epan <- flatten_plot_mean(out.epan)

  expect_gt(length(mean.gaussian), 0L)
  expect_equal(length(mean.gaussian), length(mean.epan))
  expect_gt(max(abs(mean.gaussian - mean.epan)), 1e-6)
})

test_that("condbandwidth plot helper uses method-defining x/y kernel options", {
  set.seed(9315)
  n <- 64L
  xdat <- data.frame(x = runif(n))
  ydat <- data.frame(y = runif(n))

  bw.gaussian <- npcdistbw(
    xdat = xdat,
    ydat = ydat,
    bws = c(0.22, 0.22),
    bandwidth.compute = FALSE,
    bwtype = "fixed",
    cxkertype = "gaussian",
    cykertype = "gaussian"
  )

  bw.epan <- npcdistbw(
    xdat = xdat,
    ydat = ydat,
    bws = c(0.22, 0.22),
    bandwidth.compute = FALSE,
    bwtype = "fixed",
    cxkertype = "epanechnikov",
    cykertype = "epanechnikov",
    cxkerorder = 2L,
    cykerorder = 2L
  )

  set.seed(9316)
  out.gaussian <- run_bootstrap_plot(
    bw.gaussian,
    xdat = xdat,
    ydat = ydat,
    boot_num = 9L,
    boot_method = "inid"
  )
  set.seed(9316)
  out.epan <- run_bootstrap_plot(
    bw.epan,
    xdat = xdat,
    ydat = ydat,
    boot_num = 9L,
    boot_method = "inid"
  )

  mean.gaussian <- flatten_plot_mean(out.gaussian)
  mean.epan <- flatten_plot_mean(out.epan)

  expect_gt(length(mean.gaussian), 0L)
  expect_equal(length(mean.gaussian), length(mean.epan))
  expect_gt(max(abs(mean.gaussian - mean.epan)), 1e-6)
})
