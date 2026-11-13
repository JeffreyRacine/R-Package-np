test_that("ordered absolute CDF bands retain literal block selections", {
  old <- options(np.messages = FALSE, np.plot.progress = FALSE)
  on.exit(options(old), add = TRUE)
  x <- data.frame(x = seq(-2, 2, length.out = 31),
    o = ordered(rep(c(-2, 0, 3), length.out = 31), levels = c(-2, 0, 3, 5)))
  bw <- npudistbw(dat = x, bws = c(.4, .2), bandwidth.compute = FALSE)
  for (law in c("fixed", "geom")) {
    cap <- capture_plot_provider(suppressWarnings(plot(bw, xdat = x, errors = "bootstrap",
      bootstrap = law, boot.control = np_boot_control(blocklen = 3L),
      center = "bias-corrected", band = "all", B = 49L,
      neval = 7L, output = "data", random.seed = 1312277931)))
    z <- cap$provider$d1
    flat <- plot_target(cap$result, band = "pointwise")
    expect_identical(z$inference$settings$law, law)
    expect_identical(z$inference$settings$selection$length, 3L)
    expect_false(z$inference$settings$selection$automatic)
    iv <- z$inference$bands
    expect_true(all(is.finite(unlist(iv))))
    expect_equal(flat$lower, z$bias.corrected + z$derr[, 1L], tolerance = 1e-12)
    expect_equal(flat$upper, z$bias.corrected + z$derr[, 2L], tolerance = 1e-12)
    expect_true(all(iv$pointwise[, 1L, 1L] >= iv$simultaneous[, 1L, 1L]))
    expect_true(all(iv$simultaneous[, 1L, 1L] >= iv$bonferroni[, 1L, 1L]))
    expect_true(all(iv$pointwise[, 2L, 1L] <= iv$simultaneous[, 2L, 1L]))
    expect_true(all(iv$simultaneous[, 2L, 1L] <= iv$bonferroni[, 2L, 1L]))
  }
})
