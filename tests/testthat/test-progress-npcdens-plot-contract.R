test_that("conditional bootstrap target labels format compactly", {
  fmt <- getFromNamespace(".np_plot_conditional_bootstrap_target_label", "np")

  bws <- list(
    xndim = 2L,
    yndim = 1L,
    xnames = c("x1", "x2"),
    ynames = "y"
  )

  expect_identical(fmt(bws = bws, slice.index = 0L, gradients = FALSE), "surf 1/3")
  expect_identical(fmt(bws = bws, slice.index = 1L, gradients = FALSE), "x1 1/3")
  expect_identical(fmt(bws = bws, slice.index = 3L, gradients = FALSE), "y 3/3")
  expect_identical(fmt(bws = bws, slice.index = 2L, gradients = TRUE, gradient.index = 1L), "grad x1 on x2 2/3")
})

test_that("conditional helper labels carry target context for block bootstrap phases", {
  fn <- getFromNamespace("compute.bootstrap.errors.conbandwidth", "np")

  captured <- new.env(parent = emptyenv())
  captured$prep <- character()
  captured$progress <- character()
  captured$interval <- character()

  orig_activity_begin <- getFromNamespace(".np_plot_activity_begin", "np")
  orig_activity_end <- getFromNamespace(".np_plot_activity_end", "np")
  orig_grad <- getFromNamespace(".np_inid_boot_from_conditional_gradient_local", "np")
  orig_interval <- getFromNamespace(".np_plot_bootstrap_interval_summary", "np")

  assignInNamespace(".np_plot_activity_begin", function(label) {
    captured$prep <- c(captured$prep, label)
    list(label = label)
  }, ns = "np")
  assignInNamespace(".np_plot_activity_end", function(state) invisible(NULL), ns = "np")
  assignInNamespace(".np_inid_boot_from_conditional_gradient_local", function(..., progress.label = NULL) {
    captured$progress <- c(captured$progress, progress.label)
    list(t = matrix(0, nrow = 2L, ncol = 3L), t0 = c(0, 0, 0))
  }, ns = "np")
  assignInNamespace(".np_plot_bootstrap_interval_summary", function(boot.t, t0, alpha, band.type, progress.label = NULL) {
    captured$interval <- c(captured$interval, progress.label)
    list(err = matrix(0, nrow = ncol(boot.t), ncol = 2L), all.err = list())
  }, ns = "np")
  on.exit({
    assignInNamespace(".np_plot_activity_begin", orig_activity_begin, ns = "np")
    assignInNamespace(".np_plot_activity_end", orig_activity_end, ns = "np")
    assignInNamespace(".np_inid_boot_from_conditional_gradient_local", orig_grad, ns = "np")
    assignInNamespace(".np_plot_bootstrap_interval_summary", orig_interval, ns = "np")
  }, add = TRUE)

  xdat <- data.frame(x1 = c(0, 1, 2), x2 = c(1, 2, 3))
  ydat <- data.frame(y = c(0.1, 0.2, 0.3))
  bws <- list(
    type = "fixed",
    xndim = 2L,
    yndim = 1L,
    xnames = c("x1", "x2"),
    ynames = "y",
    xdati = list(),
    ydati = list()
  )

  out <- fn(
    xdat = xdat,
    ydat = ydat,
    exdat = xdat,
    eydat = ydat,
    cdf = FALSE,
    quantreg = FALSE,
    tau = NULL,
    gradients = TRUE,
    gradient.index = 1L,
    slice.index = 2L,
    plot.errors.boot.method = "geom",
    plot.errors.boot.nonfixed = "exact",
    plot.errors.boot.blocklen = 2L,
    plot.errors.boot.num = 2L,
    plot.errors.center = "estimate",
    plot.errors.type = "all",
    plot.errors.alpha = 0.05,
    progress.target = "grad x1 on x2 2/3",
    bws = bws
  )

  expect_true(is.list(out))
  expect_true(any(grepl("Preparing plot bootstrap geom \\(grad x1 on x2 2/3\\)", captured$prep)))
  expect_true(any(grepl("Plot bootstrap \\(grad x1 on x2 2/3\\)", captured$progress)))
  expect_true(any(grepl("Constructing bootstrap all bands \\(grad x1 on x2 2/3\\)", captured$interval)))
})
