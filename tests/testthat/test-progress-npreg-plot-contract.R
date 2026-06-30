test_that("regression bootstrap target labels format compactly", {
  fmt <- getFromNamespace(".np_plot_regression_bootstrap_target_label", "np")

  bws <- list(
    ndim = 2L,
    xnames = c("x1", "x2")
  )

  expect_identical(fmt(bws = bws, slice.index = 0L, gradients = FALSE), "surf 1/2")
  expect_identical(fmt(bws = bws, slice.index = 1L, gradients = FALSE), "x1 1/2")
  expect_identical(fmt(bws = bws, slice.index = 2L, gradients = TRUE), "grad x2 2/2")
})

test_that("regression helper labels carry target context for block bootstrap phases", {
  fn <- getFromNamespace("compute.bootstrap.errors.rbandwidth", "np")

  captured <- new.env(parent = emptyenv())
  captured$prep <- character()
  captured$progress <- character()
  captured$interval <- character()

  orig_activity_begin <- getFromNamespace(".np_plot_activity_begin", "np")
  orig_activity_end <- getFromNamespace(".np_plot_activity_end", "np")
  orig_reg_boot <- getFromNamespace(".np_inid_boot_from_regression", "np")
  orig_interval <- getFromNamespace(".np_plot_bootstrap_interval_summary", "np")

  assignInNamespace(".np_plot_activity_begin", function(label) {
    captured$prep <- c(captured$prep, label)
    list(label = label)
  }, ns = "np")
  assignInNamespace(".np_plot_activity_end", function(state) invisible(NULL), ns = "np")
  assignInNamespace(".np_inid_boot_from_regression", function(..., prep.label = NULL, progress.label = NULL) {
    captured$progress <- c(captured$progress, prep.label, progress.label)
    list(t = matrix(0, nrow = 2L, ncol = 3L), t0 = c(0, 0, 0))
  }, ns = "np")
  assignInNamespace(".np_plot_bootstrap_interval_summary", function(boot.t, t0, alpha, band.type, progress.label = NULL) {
    captured$interval <- c(captured$interval, progress.label)
    list(err = matrix(0, nrow = ncol(boot.t), ncol = 2L), all.err = list())
  }, ns = "np")
  on.exit({
    assignInNamespace(".np_plot_activity_begin", orig_activity_begin, ns = "np")
    assignInNamespace(".np_plot_activity_end", orig_activity_end, ns = "np")
    assignInNamespace(".np_inid_boot_from_regression", orig_reg_boot, ns = "np")
    assignInNamespace(".np_plot_bootstrap_interval_summary", orig_interval, ns = "np")
  }, add = TRUE)

  xdat <- data.frame(x1 = c(0, 1, 2), x2 = c(1, 2, 3))
  bws <- list(
    type = "fixed",
    ndim = 2L,
    xnames = c("x1", "x2"),
    xdati = list(
      icon = c(TRUE, TRUE),
      iord = c(FALSE, FALSE),
      iuno = c(FALSE, FALSE)
    )
  )

  out <- fn(
    xdat = xdat,
    ydat = c(1, 2, 3),
    exdat = xdat,
    gradients = TRUE,
    gradient.order = 1L,
    slice.index = 2L,
    plot.errors.boot.method = "geom",
    plot.errors.boot.blocklen = 2L,
    plot.errors.boot.num = 2L,
    plot.errors.center = "estimate",
    plot.errors.type = "all",
    plot.errors.alpha = 0.05,
    progress.target = "grad x2 2/2",
    bws = bws
  )

  expect_true(is.list(out))
  expect_true(any(grepl("Preparing plot bootstrap geom \\(grad x2 2/2\\)", captured$prep, fixed = FALSE)))
  expect_true(any(grepl("Preparing plot bootstrap geom \\(grad x2 2/2\\)", captured$progress, fixed = FALSE)))
  expect_true(any(grepl("Plot bootstrap \\(grad x2 2/2\\)", captured$progress, fixed = FALSE)))
  expect_true(any(grepl("Constructing bootstrap all bands \\(grad x2 2/2\\)", captured$interval, fixed = FALSE)))
})

test_that("regression helper labels carry target context for wild bootstrap phases", {
  fn <- getFromNamespace("compute.bootstrap.errors.rbandwidth", "np")

  captured <- new.env(parent = emptyenv())
  captured$prep <- character()
  captured$progress <- character()
  captured$interval <- character()

  orig_activity_begin <- getFromNamespace(".np_plot_activity_begin", "np")
  orig_activity_end <- getFromNamespace(".np_plot_activity_end", "np")
  orig_wild <- getFromNamespace(".np_plot_boot_from_hat_wild", "np")
  orig_interval <- getFromNamespace(".np_plot_bootstrap_interval_summary", "np")
  orig_npreghat <- getFromNamespace("npreghat", "np")

  assignInNamespace(".np_plot_activity_begin", function(label) {
    captured$prep <- c(captured$prep, label)
    list(label = label)
  }, ns = "np")
  assignInNamespace(".np_plot_activity_end", function(state) invisible(NULL), ns = "np")
  assignInNamespace(".np_plot_boot_from_hat_wild", function(H, ydat, fit.mean, B, wild, progress.label = NULL) {
    captured$progress <- c(captured$progress, progress.label)
    list(t = matrix(0, nrow = 2L, ncol = 3L), t0 = c(0, 0, 0))
  }, ns = "np")
  assignInNamespace(".np_plot_bootstrap_interval_summary", function(boot.t, t0, alpha, band.type, progress.label = NULL) {
    captured$interval <- c(captured$interval, progress.label)
    list(err = matrix(0, nrow = ncol(boot.t), ncol = 2L), all.err = list())
  }, ns = "np")
  assignInNamespace("npreghat", function(...) {
    args <- list(...)
    if (identical(args$output, "apply")) {
      return(rep(0, nrow(args$txdat)))
    }
    matrix(0, nrow = nrow(args$exdat), ncol = nrow(args$txdat))
  }, ns = "np")
  on.exit({
    assignInNamespace(".np_plot_activity_begin", orig_activity_begin, ns = "np")
    assignInNamespace(".np_plot_activity_end", orig_activity_end, ns = "np")
    assignInNamespace(".np_plot_boot_from_hat_wild", orig_wild, ns = "np")
    assignInNamespace(".np_plot_bootstrap_interval_summary", orig_interval, ns = "np")
    assignInNamespace("npreghat", orig_npreghat, ns = "np")
  }, add = TRUE)

  xdat <- data.frame(x1 = c(0, 1, 2), x2 = c(1, 2, 3))
  bws <- list(
    type = "fixed",
    ndim = 2L,
    xnames = c("x1", "x2"),
    xdati = list(
      icon = c(TRUE, TRUE),
      iord = c(FALSE, FALSE),
      iuno = c(FALSE, FALSE)
    )
  )

  out <- fn(
    xdat = xdat,
    ydat = c(1, 2, 3),
    exdat = xdat,
    gradients = TRUE,
    gradient.order = 1L,
    slice.index = 1L,
    plot.errors.boot.method = "wild",
    plot.errors.boot.wild = "rademacher",
    plot.errors.boot.blocklen = NULL,
    plot.errors.boot.num = 2L,
    plot.errors.center = "estimate",
    plot.errors.type = "all",
    plot.errors.alpha = 0.05,
    progress.target = "grad x1 1/2",
    bws = bws
  )

  expect_true(is.list(out))
  expect_true(any(grepl("Preparing plot bootstrap wild \\(grad x1 1/2\\)", captured$prep, fixed = FALSE)))
  expect_true(any(grepl("Plot bootstrap \\(grad x1 1/2\\)", captured$progress, fixed = FALSE)))
  expect_true(any(grepl("Constructing bootstrap all bands \\(grad x1 1/2\\)", captured$interval, fixed = FALSE)))
})
