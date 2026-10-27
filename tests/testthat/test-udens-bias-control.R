test_that("density controls are explicit and do not change regression defaults", {
  control <- np_udens_bias_control()
  expect_s3_class(control, "np_udens_bias_control")
  expect_identical(control$construction, "positive-covariance")
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
  out <- plot(bw, xdat = x, errors = "bootstrap", bootstrap = "inid",
    center = "bias-corrected", band = "all", B = 49L, neval = 1L,
    output = "data", bias.control = np_udens_bias_control())[[1L]]
  expect_identical(.Random.seed, rng)
  expect_identical(out$inference$settings$primary, "absolute")
  expect_identical(out$inference$bands$pointwise, out$inference$bands$simultaneous)
  expect_identical(out$inference$bands$pointwise, out$inference$bands$bonferroni)
  expect_error(plot(bw, xdat = x, errors = "bootstrap", bootstrap = "fixed",
    center = "bias-corrected", band = "all", B = 49L, neval = 1L,
    boot.control = np_boot_control(blocklen = 2L), output = "data",
    bias.control = np_udens_bias_control()), "bootstrap='inid'", fixed = TRUE)
})
