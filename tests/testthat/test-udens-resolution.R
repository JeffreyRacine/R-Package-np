test_that("density low-draw warnings aggregate across slices", {
  old <- options(np.messages = FALSE, np.plot.progress = FALSE)
  on.exit(options(old), add = TRUE)
  x <- data.frame(x = sin(1:20), y = cos(1:20))
  bw <- npudensbw(dat = x, bws = c(.6, .6), bandwidth.compute = FALSE)
  warnings <- list()
  set.seed(716)
  seed <- .Random.seed
  out <- withCallingHandlers(plot(bw, xdat = x, errors = "bootstrap",
    center = "bias-corrected", band = "all", B = 1L, neval = 3L,
    perspective = FALSE, output = "data"), warning = function(w) {
      warnings[[length(warnings) + 1L]] <<- w
      invokeRestart("muffleWarning")
    })
  expect_identical(.Random.seed, seed)
  expect_length(warnings, 1L)
  expect_s3_class(warnings[[1L]], "np_band_resolution_warning")
  for (z in out) {
    expect_named(z$inference$bands, c("pointwise", "simultaneous", "bonferroni"))
    expect_true(all(is.finite(unlist(z$inference$bands))))
    expect_true(all(z$inference$resolution$resolution_limited))
    expect_true(any(z$inference$resolution$band == "scale-simultaneous"))
  }
})
