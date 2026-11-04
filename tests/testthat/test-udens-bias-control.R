test_that("density controls are explicit and do not change regression defaults", {
  control <- np_udens_bias_control()
  expect_s3_class(control, "np_udens_bias_control")
  expect_identical(control$construction, "joint-scale")
  expect_identical(control$version, 3L)
  expect_identical(np_udens_bias_control("positive-covariance")$version, 1L)
  expect_identical(control$calibration, "absolute")
  expect_identical(np_bias_control()$iterations, 16L)
  control$version <- 9L
  expect_error(npRmpi:::.np_udens_plot_prepare(control, NULL, NULL,
    "bias-corrected", "bootstrap", "inid", "all"), "np_udens_bias_control")
  expect_null(npRmpi:::.np_udens_plot_prepare(NULL, NULL, NULL,
    "estimate", "asymptotic", "inid", "pointwise"))
})

test_that("density correction retains literal law and singleton contracts", {
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  set.seed(731120001)
  x <- data.frame(x = rnorm(40))
  bw <- npudensbw(dat = x, bws = .5, bandwidth.compute = FALSE)
  rng <- .Random.seed
  grDevices::pdf(file = NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  out <- plot(bw, xdat = x, errors = "bootstrap", bootstrap = "inid",
    center = "bias-corrected", band = "all", B = 49L, neval = 1L,
    output = "plot-data", bias_control = np_udens_bias_control())
  expect_identical(.Random.seed, rng)
  expect_identical(out$inference[[1L]]$primary, "absolute")
  expect_identical(plot_band(out, "pointwise"), plot_band(out, "simultaneous"))
  expect_identical(plot_band(out, "pointwise"), plot_band(out, "bonferroni"))
  for (law in c("fixed", "geom")) {
    block <- plot(bw, xdat = x, errors = "bootstrap", bootstrap = law,
      center = "bias-corrected", band = "all", B = 49L, neval = 1L,
      boot_control = np_boot_control(blocklen = 2L), output = "data",
      bias_control = np_udens_bias_control())
    expect_identical(block$inference[[1L]]$resampling_law, law)
    expect_true(all(is.finite(as.matrix(block$panels[[1L]]$intervals[c("lower", "upper")]))))
    expect_identical(plot_band(block, "pointwise"), plot_band(block, "simultaneous"))
    expect_identical(plot_band(block, "pointwise"), plot_band(block, "bonferroni"))
    expect_identical(.Random.seed, rng)
  }
})


test_that("density construction names retain saved-control meanings", {
  old <- options(np.messages = FALSE, np.plot.progress = FALSE)
  on.exit(options(old), add = TRUE)
  x <- data.frame(x = seq(-2, 2, length.out = 12L))
  bw <- npudensbw(dat = x, bws = .6, bandwidth.compute = FALSE)
  prepare <- getFromNamespace(".np_udens_plot_prepare", "npRmpi")
  for (version in 1:3) {
    saved <- structure(list(construction = "positive-covariance",
      calibration = "absolute", version = as.integer(version)),
      class = "np_udens_bias_control")
    state <- prepare(saved, bw, x, "bias-corrected", "bootstrap", "inid", "pointwise")
    expect_identical(state$plot.control.version, as.integer(version))
  }
  control <- np_udens_bias_control()
  control$version <- 1L
  expect_error(prepare(control, bw, x, "bias-corrected", "bootstrap", "inid", "pointwise"),
    "joint-scale density controls require version 3", fixed = TRUE)
  expect_null(prepare(np_udens_bias_control("legacy"), bw, x,
    "bias-corrected", "bootstrap", "inid", "pointwise"))
})
