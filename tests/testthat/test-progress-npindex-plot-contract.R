test_that("single-index bootstrap target labels format compactly", {
  fmt <- getFromNamespace(".np_plot_singleindex_bootstrap_target_label", "np")

  expect_null(fmt(gradients = FALSE))
  expect_null(fmt(gradients = TRUE))
})

test_that("single-index helper labels carry target context for bootstrap phases", {
  fn <- getFromNamespace("compute.bootstrap.errors.sibandwidth", "np")

  captured <- new.env(parent = emptyenv())
  captured$prep <- character()
  captured$progress <- character()
  captured$interval <- character()

  orig_activity_begin <- getFromNamespace(".np_plot_activity_begin", "np")
  orig_activity_end <- getFromNamespace(".np_plot_activity_end", "np")
  orig_index <- getFromNamespace(".np_inid_boot_from_index", "np")
  orig_interval <- getFromNamespace(".np_plot_bootstrap_interval_summary", "np")

  assignInNamespace(".np_plot_activity_begin", function(label) {
    captured$prep <- c(captured$prep, label)
    list(label = label)
  }, ns = "np")
  assignInNamespace(".np_plot_activity_end", function(state) invisible(NULL), ns = "np")
  assignInNamespace(".np_inid_boot_from_index", function(..., progress.label = NULL) {
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
    assignInNamespace(".np_inid_boot_from_index", orig_index, ns = "np")
    assignInNamespace(".np_plot_bootstrap_interval_summary", orig_interval, ns = "np")
  }, add = TRUE)

  xdat <- data.frame(x1 = c(0, 1, 2), x2 = c(1, 2, 3))
  ydat <- c(0.1, 0.2, 0.3)
  bws <- list(
    type = "fixed",
    beta = c(1, 1)
  )

  out <- fn(
    xdat = xdat,
    ydat = ydat,
    gradients = FALSE,
    plot.errors.boot.method = "geom",
    plot.errors.boot.blocklen = 2L,
    plot.errors.boot.num = 2L,
    plot.errors.center = "estimate",
    plot.errors.type = "all",
    plot.errors.alpha = 0.05,
    bws = bws
  )

  expect_true(is.list(out))
  expect_true(any(grepl("^Preparing plot bootstrap geom$", captured$prep)))
  expect_true(any(grepl("^Plot bootstrap$", captured$progress)))
  expect_true(any(grepl("^Constructing bootstrap all bands$", captured$interval)))
})
