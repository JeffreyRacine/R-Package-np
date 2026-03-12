test_that("bootstrap interval summary avoids recursive pointwise recomputation for all-type requests", {
  skip_if_not_installed("npRmpi")

  interval_summary <- getFromNamespace(".np_plot_bootstrap_interval_summary", "npRmpi")
  quantile_bounds <- getFromNamespace("compute.bootstrap.quantile.bounds", "npRmpi")

  set.seed(20260310L)
  boot.t <- matrix(rnorm(128L * 5L), nrow = 128L, ncol = 5L)
  t0 <- colMeans(boot.t) + seq_len(ncol(boot.t)) / 10
  alpha <- 0.1

  assign("np_plot_band_type_calls_test", character(), envir = .GlobalEnv)
  trace(
    "compute.bootstrap.quantile.bounds",
    where = asNamespace("npRmpi"),
    print = FALSE,
    tracer = quote({
      assign(
        "np_plot_band_type_calls_test",
        c(get("np_plot_band_type_calls_test", envir = .GlobalEnv, inherits = FALSE), band.type),
        envir = .GlobalEnv
      )
    })
  )
  on.exit({
    untrace("compute.bootstrap.quantile.bounds", where = asNamespace("npRmpi"))
    rm("np_plot_band_type_calls_test", envir = .GlobalEnv)
  }, add = TRUE)

  out <- interval_summary(
    boot.t = boot.t,
    t0 = t0,
    alpha = alpha,
    band.type = "all"
  )

  expect_identical(
    get("np_plot_band_type_calls_test", envir = .GlobalEnv, inherits = FALSE),
    c("all", "simultaneous")
  )

  pointwise <- quantile_bounds(boot.t = boot.t, alpha = alpha, band.type = "pointwise", warn.coverage = FALSE)
  bonferroni <- quantile_bounds(boot.t = boot.t, alpha = alpha, band.type = "bonferroni", warn.coverage = FALSE)
  simultaneous <- quantile_bounds(boot.t = boot.t, alpha = alpha, band.type = "simultaneous", warn.coverage = FALSE)

  expect_equal(out$err, cbind(t0 - pointwise[, 1L], pointwise[, 2L] - t0))
  expect_equal(out$all.err$pointwise, cbind(t0 - pointwise[, 1L], pointwise[, 2L] - t0))
  expect_equal(out$all.err$bonferroni, cbind(t0 - bonferroni[, 1L], bonferroni[, 2L] - t0))
  expect_equal(out$all.err$simultaneous, cbind(t0 - simultaneous[, 1L], simultaneous[, 2L] - t0))
})

test_that("pmzsd helper matches covariance-diagonal standard deviations", {
  skip_if_not_installed("npRmpi")

  boot_col_sds <- getFromNamespace(".np_plot_bootstrap_col_sds", "npRmpi")

  set.seed(20260311L)
  boot.t <- matrix(rnorm(257L * 11L), nrow = 257L, ncol = 11L)

  expect_equal(
    boot_col_sds(boot.t),
    sqrt(diag(cov(boot.t))),
    tolerance = 1e-12
  )
})

test_that("sibandwidth bootstrap route uses the shared interval summary helper", {
  skip_if_not_installed("npRmpi")

  fn <- getFromNamespace("compute.bootstrap.errors.sibandwidth", "npRmpi")
  txt <- paste(deparse(body(fn), width.cutoff = 500L), collapse = "\n")

  expect_match(txt, "\\.np_plot_bootstrap_interval_summary\\(", perl = TRUE)
  expect_false(
    grepl("boot\\.bounds <- compute\\.bootstrap\\.quantile\\.bounds\\(", txt, perl = TRUE),
    info = "sibandwidth should not maintain a private pointwise/all interval path"
  )
  expect_false(
    grepl("boot\\.all\\.bounds <- compute\\.bootstrap\\.quantile\\.bounds\\(", txt, perl = TRUE),
    info = "sibandwidth should reuse the shared all-band summary path"
  )
})
