library(np)

with_npreg_plot_device <- function(expr) {
  tf <- tempfile(fileext = ".pdf")
  grDevices::pdf(tf)
  on.exit({
    grDevices::dev.off()
    unlink(tf)
  }, add = TRUE)
  force(expr)
}

test_that("npreg adaptive-nn wild gradient plot-data stays available in serial", {
  set.seed(42)
  n <- 120L
  x <- sort(stats::runif(n))
  y <- sin(2 * pi * x) + stats::rnorm(n, sd = 0.1)
  xdat <- data.frame(x = x)

  bw <- npregbw(
    xdat = xdat,
    ydat = y,
    regtype = "ll",
    bwtype = "adaptive_nn",
    bws = 35L,
    bandwidth.compute = FALSE
  )
  fit <- npreg(
    bws = bw,
    txdat = xdat,
    tydat = y,
    warn.glp.gradient = FALSE
  )

  out <- with_npreg_plot_device(suppressWarnings(plot(
    fit,
    xdat = xdat,
    ydat = y,
    plot.behavior = "data",
    perspective = FALSE,
    neval = 40L,
    plot.errors.method = "bootstrap",
    plot.errors.boot.method = "wild",
    plot.errors.boot.num = 39L,
    plot.errors.type = "all",
    gradients = TRUE
  )))

  expect_type(out, "list")
  expect_s3_class(out[[1L]], "npregression")
  expect_length(out[[1L]]$mean, 40L)
})
