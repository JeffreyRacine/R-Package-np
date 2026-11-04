test_that("explicit joint-scale density surfaces expose signed endpoint offsets", {
  old <- options(np.messages = FALSE, np.plot.progress = FALSE)
  on.exit(options(old), add = TRUE)
  x <- data.frame(x1 = c(-2, -.5, .2, 1.4, 2.1),
                  x2 = c(1.5, -1.2, .4, 2.3, -.8))
  bw <- npudensbw(dat = x, bws = c(.8, .9), bandwidth.compute = FALSE)
  fit <- npudens(bws = bw, tdat = x)
  set.seed(731120021)
  rng <- .Random.seed
  for (method in c("asymptotic", "bootstrap")) {
    args <- list(x = fit, errors = method, center = "bias-corrected",
                 band = "all", neval = 2L, output = "data",
                 bias.control = np_udens_bias_control())
    args[[if (method == "asymptotic") "gaussian.draws" else "B"]] <- 199L
    if (method == "asymptotic")
      expect_warning(expect_warning(out <- do.call(plot, args),
        "gaussian.draws=199", fixed = TRUE), "gaussian.draws=199", fixed = TRUE)
    else out <- do.call(plot, args)
    r <- out$d1
    expect_identical(.Random.seed, rng)
    expect_identical(r$inference$settings$law,
                     if (method == "asymptotic") "gaussian" else "inid")
    expect_equal(range(r$eval[[1L]]), range(x[[1L]]), tolerance = 0)
    expect_equal(range(r$eval[[2L]]), range(x[[2L]]), tolerance = 0)
    expect_equal(r$bias.corrected + r$derr[, 1L],
                 r$inference$bands$pointwise[, 1L, 1L], tolerance = 1e-14)
    expect_equal(r$bias.corrected + r$derr[, 2L],
                 r$inference$bands$pointwise[, 2L, 1L], tolerance = 1e-14)
  }
})

test_that("underresolved corrected bootstrap Bonferroni is explicit", {
  old <- options(np.messages = FALSE, np.plot.progress = FALSE)
  on.exit(options(old), add = TRUE)
  x <- data.frame(x = seq(-2, 2, length.out = 10))
  bw <- npudensbw(dat = x, bws = .5, bandwidth.compute = FALSE)
  set.seed(731120022)
  rng <- .Random.seed
  expect_warning(out <- plot(bw, xdat = x, errors = "bootstrap",
    center = "bias-corrected", band = "bonferroni", B = 49L, neval = 3L,
    output = "data"), "limited draw resolution", fixed = TRUE)
  expect_true(all(is.finite(out$d1$inference$bands$bonferroni)))
  expect_true(any(out$d1$inference$resolution$resolution_limited))
  expect_identical(.Random.seed, rng)
  # The same explicit law and draw count still admit the distinct maximum band.
  expect_s3_class(plot(bw, xdat = x, errors = "bootstrap",
    center = "bias-corrected", band = "simultaneous", B = 49L, neval = 3L,
    output = "data")$d1, "npdensity")
})
