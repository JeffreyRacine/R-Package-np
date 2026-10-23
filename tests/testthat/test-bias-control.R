test_that("bias controls preserve explicit construction and auxiliary basis", {
  control <- np_bias_control()
  expect_s3_class(control, "np_bias_control")
  expect_identical(control$iterations, 16L)
  expect_true(control$bernstein.basis)
  expect_identical(np_boot_control()$wild, "mammen")
  expect_identical(np_boot_control(wild="rademacher")$wild, "rademacher")
  expect_error(np_bias_control(degree=3, bernstein.basis=FALSE), "Bernstein|bernstein")
  expect_error(np_bias_control(iterations=1.2), "positive integer")
  expect_error(np_bias_control(p0.scale=4, p0.bandwidth=.2), "mutually exclusive")
  expect_error(np_bias_control(construction="legacy", degree=5), "do not apply")
  control$construction <- "typo"
  expect_error(npRmpi:::.np_bias_validate_control(control), "invalid np_bias_control")
})

test_that("corrected plot families nest and singleton limits coincide", {
  if (!spawn_mpi_slaves(n=1L)) skip("MPI pool unavailable in this test context")
  on.exit(close_mpi_slaves(force=TRUE), add=TRUE)
  old <- options(np.messages=FALSE)
  on.exit(options(old), add=TRUE)
  set.seed(27051)
  x <- data.frame(x=runif(80))
  y <- sin(3*x$x) + rnorm(80, sd=.2)
  bw <- npregbw(xdat=x, ydat=y, bws=.23, bandwidth.compute=FALSE,
                regtype="lp", degree=2, bernstein.basis=TRUE)
  fit <- npreg(bw, txdat=x, tydat=y)
  original <- fit
  rng <- .Random.seed
  out <- plot(fit, errors="asymptotic", center="bias-corrected", band="all",
              gaussian.draws=63L, neval=3L, output="data")[[1L]]
  expect_identical(fit, original)
  expect_identical(.Random.seed, rng)
  expect_true(out$inference$settings$bernstein.basis)
  expect_identical(out$inference$settings$main$degree, bw$degree)
  b <- out$inference$bands
  expect_true(all(b$bonferroni[,1,] <= b$simultaneous[,1,]))
  expect_true(all(b$simultaneous[,1,] <= b$pointwise[,1,]))
  expect_true(all(b$pointwise[,2,] <= b$simultaneous[,2,]))
  expect_true(all(b$simultaneous[,2,] <= b$bonferroni[,2,]))
  state <- npRmpi:::.np_bias_regression_prepare(bw, x, y)
  for (law in c("gaussian", "mammen", "rademacher", "inid", "fixed", "geom")) {
    z <- npRmpi:::.np_with_seed(42L, npRmpi:::.np_bias_inference(
      state, x[1L,,drop=FALSE], "mean", .05, 63L, law))
    expect_identical(z$bands$pointwise, z$bands$simultaneous)
    expect_identical(z$bands$pointwise, z$bands$bonferroni)
  }
})

test_that("iterated intervals describe maximum calibration without relabeling legacy bands", {
  label <- npRmpi:::.np_plot_all_band_legend_labels
  descriptor <- npRmpi:::.np_plot_variability_descriptor
  expect_match(label("bootstrap")[["simultaneous"]], "rank-based", fixed=TRUE)
  expect_match(label("bootstrap", "iterated")[["simultaneous"]], "maximum", fixed=TRUE)
  expect_match(label("bootstrap", "iterated")[["bonferroni"]], "envelope", fixed=TRUE)
  text <- descriptor("bootstrap", "simultaneous", .05, "bias-corrected",
                     construction="iterated")$candidates
  expect_true(all(grepl("maximum", text, fixed=TRUE)))
  expect_false(any(grepl("rank", text, fixed=TRUE)))
  text <- descriptor("asymptotic", "pointwise", .05, "bias-corrected",
                     construction="iterated")$candidates
  expect_true(all(grepl("bias-corrected", text, fixed=TRUE)))
})
