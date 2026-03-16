test_that("smooth coefficient bootstrap target labels format compactly", {
  fmt <- getFromNamespace(".np_plot_scoef_bootstrap_target_label", "np")

  bws <- list(
    xndim = 1L,
    zndim = 1L,
    xnames = "x",
    znames = "z"
  )

  expect_identical(fmt(bws = bws, slice.index = 0L), "surf 1/2")
  expect_identical(fmt(bws = bws, slice.index = 1L), "x 1/2")
  expect_identical(fmt(bws = bws, slice.index = 2L), "z 2/2")
})

test_that("smooth coefficient helper labels carry target context for bootstrap phases", {
  fn <- getFromNamespace("compute.bootstrap.errors.scbandwidth", "np")

  captured <- new.env(parent = emptyenv())
  captured$prep <- character()
  captured$progress <- character()
  captured$interval <- character()

  orig_activity_begin <- getFromNamespace(".np_plot_activity_begin", "np")
  orig_activity_end <- getFromNamespace(".np_plot_activity_end", "np")
  orig_scoef <- getFromNamespace(".np_inid_boot_from_scoef", "np")
  orig_interval <- getFromNamespace(".np_plot_bootstrap_interval_summary", "np")

  assignInNamespace(".np_plot_activity_begin", function(label) {
    captured$prep <- c(captured$prep, label)
    list(label = label)
  }, ns = "np")
  assignInNamespace(".np_plot_activity_end", function(state) invisible(NULL), ns = "np")
  assignInNamespace(".np_inid_boot_from_scoef", function(..., progress.label = NULL) {
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
    assignInNamespace(".np_inid_boot_from_scoef", orig_scoef, ns = "np")
    assignInNamespace(".np_plot_bootstrap_interval_summary", orig_interval, ns = "np")
  }, add = TRUE)

  xdat <- data.frame(x = c(0, 1, 2))
  zdat <- data.frame(z = c(1, 2, 3))
  ydat <- c(0.1, 0.2, 0.3)
  bws <- list(
    xndim = 1L,
    zndim = 1L,
    xnames = "x",
    znames = "z"
  )

  out <- fn(
    xdat = xdat,
    ydat = ydat,
    zdat = zdat,
    exdat = xdat,
    ezdat = zdat,
    gradients = FALSE,
    slice.index = 2L,
    progress.target = "z 2/2",
    plot.errors.boot.method = "geom",
    plot.errors.boot.blocklen = 2L,
    plot.errors.boot.num = 2L,
    plot.errors.center = "estimate",
    plot.errors.type = "all",
    plot.errors.alpha = 0.05,
    bws = bws
  )

  expect_true(is.list(out))
  expect_true(any(grepl("Preparing plot bootstrap geom \\(z 2/2\\)", captured$prep)))
  expect_true(any(grepl("Plot bootstrap \\(z 2/2\\)", captured$progress)))
  expect_true(any(grepl("Constructing bootstrap all bands \\(z 2/2\\)", captured$interval)))
})
