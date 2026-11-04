# Raw process/draw invariants are provider contracts; the public result is
# constructed and validated by capture_plot_provider in every call.
test_that("CDF raw and protected tails preserve their same-call absolute control", {
  old <- options(np.messages = FALSE, np.plot.progress = FALSE)
  on.exit(options(old), add = TRUE)
  x <- data.frame(x = seq(-2.5, 2.5, length.out = 31L) + .1*sin(1:31))
  bw <- npudistbw(dat = x, bws = .45, bandwidth.compute = FALSE)
  for (law in c("inid", "fixed", "geom")) {
    call <- list(x = bw, xdat = x, center = "bias-corrected", errors = "bootstrap",
      bootstrap = law, B = 199L, band = "all", neval = 5L, output = "data",
      random_seed = 99102905)
    block <- if (law == "inid") NULL else 3L
    call$boot_control <- np_boot_control(blocklen = block)
    absolute <- capture_plot_provider(suppressWarnings(do.call(plot, call)))$provider$d1$inference
    call$boot_control <- np_boot_control(blocklen = block, calibration = "separate-tails")
    raw <- capture_plot_provider(suppressWarnings(do.call(plot, call)))$provider$d1$inference
    call$boot_control <- np_boot_control(blocklen = block, calibration = "tail-envelope")
    protected <- capture_plot_provider(suppressWarnings(do.call(plot, call)))$provider$d1$inference
    expect_identical(raw$absolute.intervals, absolute$bands)
    expect_identical(raw$rng, absolute$rng)
    expect_identical(raw$center, absolute$center)
    expect_identical(raw$se, absolute$se)
    expect_identical(protected$raw.intervals, raw$bands)
    for (band in names(raw$bands)) {
      expect_identical(protected$bands[[band]][, 1L, ],
        pmin(raw$bands[[band]][, 1L, ], absolute$bands[[band]][, 1L, ]))
      expect_identical(protected$bands[[band]][, 2L, ],
        pmax(raw$bands[[band]][, 2L, ], absolute$bands[[band]][, 2L, ]))
    }
  }
})

test_that("CDF side-specific score protection does not impose symmetry", {
  center <- c(.02, .15, .6)
  widths <- list(pointwise = .001, simultaneous = .002, bonferroni = .003)
  intervals <- lapply(widths, function(w) {
    z <- array(0, c(3L, 2L, 1L)); z[, 1L, 1L] <- center-w
    z[, 2L, 1L] <- center+w; z
  })
  z <- .np_udist_side_envelope(intervals, center, 31L, 1L, .05, TRUE)$intervals
  expect_true(any(abs((center-z$simultaneous[, 1L, 1L]) -
                       (z$simultaneous[, 2L, 1L]-center)) > 1e-4))
  expect_true(all(z$simultaneous[, 1L, ] <= z$pointwise[, 1L, ]))
  expect_true(all(z$simultaneous[, 2L, ] >= z$pointwise[, 2L, ]))
  expect_true(all(z$bonferroni[, 1L, ] <= z$simultaneous[, 1L, ]))
  expect_true(all(z$bonferroni[, 2L, ] >= z$simultaneous[, 2L, ]))
  expect_identical(np_udist_bias_control(), structure(
    list(construction = "joint-process", version = 1L), class = "np_udist_bias_control"))
  expect_identical(np_udist_bias_control(calibration = "score-envelope")$version, 2L)
  expect_error(np_udist_bias_control("legacy", calibration = "score-envelope"), "joint-process")
})

test_that("new CDF score controls are rejected at other public consumers", {
  control <- np_udist_bias_control(calibration = "score-envelope")
  expect_error(plot(structure(list(), class = "npcopula"), errors = "bootstrap", B = 19L, output = "data",
                    bias_control = control), "requires supported distribution inference")
  expect_error(.np_udist_score_scope(control, FALSE), "requires supported distribution inference")
  expect_invisible(.np_udist_score_scope(np_udist_bias_control(), FALSE))
  expect_invisible(.np_udist_score_scope(control, TRUE))
})
