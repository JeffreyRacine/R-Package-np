test_that("distribution correction uses the public center and signed surface offsets", {
  old <- options(np.messages = FALSE, np.plot.progress = FALSE)
  on.exit(options(old), add = TRUE)
  x <- data.frame(x1 = seq(-2, 2, length.out = 17), x2 = cos(seq_len(17)))
  bw <- npudistbw(dat = x, bws = c(.6, .5), bandwidth.compute = FALSE)
  for (method in c("asymptotic", "bootstrap")) {
    args <- list(x = bw, xdat = x, errors = method, center = "bias-corrected",
      band = "all", neval = 7L, output = "data")
    args[[if (method == "asymptotic") "gaussian.draws" else "B"]] <- 49L
    z <- suppressWarnings(do.call(plot, args))$d1
    expect_identical(z$inference$settings$construction, "cdf-absolute-process-v1")
    expect_identical(z$inference$settings$critical.units, "CDF")
    expect_equal(nrow(z$eval), 49L)
    expect_equal(range(z$eval[[1L]]), range(x$x1), tolerance = 0)
    expect_equal(range(z$eval[[2L]]), range(x$x2), tolerance = 0)
    expect_equal(z$bias.corrected + z$derr[,1L], z$inference$bands$pointwise[,1L,1L], tolerance = 1e-12)
    expect_equal(z$bias.corrected + z$derr[,2L], z$inference$bands$pointwise[,2L,1L], tolerance = 1e-12)
    args$neval <- 1L
    one <- suppressWarnings(do.call(plot, args))$d1
    expect_identical(one$inference$bands$pointwise, one$inference$bands$simultaneous)
    expect_identical(one$inference$bands$pointwise, one$inference$bands$bonferroni)
  }
})

test_that("distribution controls reject altered objects and retain regression settings", {
  prepare <- getFromNamespace(".np_udist_plot_prepare", "np")
  expect_identical(np_udist_bias_control()$version, 1L)
  expect_identical(np_udist_bias_control()$construction, "joint-process")
  expect_identical(np_bias_control()$iterations, 16L)
  expect_null(prepare(NULL, NULL, NULL, "estimate", "asymptotic", "inid", "pointwise"))
  expect_null(prepare(np_udist_bias_control("legacy"), NULL, NULL,
    "bias-corrected", "bootstrap", "inid", "pointwise"))
  expect_error(prepare(np_udist_bias_control("legacy"), NULL, NULL,
    "bias-corrected", "asymptotic", "inid", "pointwise"), "legacy distribution")
  bad <- np_udist_bias_control(); bad$version <- 2L
  expect_error(prepare(bad, NULL, NULL, "bias-corrected", "bootstrap", "inid", "pointwise"), "must be created")
})
