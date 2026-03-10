test_that("bootstrap interval summary reuses pointwise bounds for all-type requests", {
  skip_if_not_installed("np")

  interval_summary <- getFromNamespace(".np_plot_bootstrap_interval_summary", "np")
  quantile_bounds <- getFromNamespace("compute.bootstrap.quantile.bounds", "np")

  set.seed(20260310L)
  boot.t <- matrix(rnorm(128L * 5L), nrow = 128L, ncol = 5L)
  t0 <- colMeans(boot.t) + seq_len(ncol(boot.t)) / 10
  alpha <- 0.1

  pointwise_calls <- 0L
  assign("np_plot_pointwise_calls_test", pointwise_calls, envir = .GlobalEnv)
  trace(
    "compute.bootstrap.quantile.bounds",
    where = asNamespace("np"),
    print = FALSE,
    tracer = quote({
      if (identical(band.type, "pointwise"))
        assign(
          "np_plot_pointwise_calls_test",
          get("np_plot_pointwise_calls_test", envir = .GlobalEnv, inherits = FALSE) + 1L,
          envir = .GlobalEnv
        )
    })
  )
  on.exit({
    untrace("compute.bootstrap.quantile.bounds", where = asNamespace("np"))
    rm("np_plot_pointwise_calls_test", envir = .GlobalEnv)
  }, add = TRUE)

  out <- interval_summary(
    boot.t = boot.t,
    t0 = t0,
    alpha = alpha,
    band.type = "all"
  )

  expect_identical(get("np_plot_pointwise_calls_test", envir = .GlobalEnv, inherits = FALSE), 1L)

  pointwise <- quantile_bounds(boot.t = boot.t, alpha = alpha, band.type = "pointwise", warn.coverage = FALSE)
  bonferroni <- quantile_bounds(boot.t = boot.t, alpha = alpha, band.type = "bonferroni", warn.coverage = FALSE)
  simultaneous <- quantile_bounds(boot.t = boot.t, alpha = alpha, band.type = "simultaneous", warn.coverage = FALSE)

  expect_equal(out$err, cbind(t0 - pointwise[, 1L], pointwise[, 2L] - t0))
  expect_equal(out$all.err$pointwise, cbind(t0 - pointwise[, 1L], pointwise[, 2L] - t0))
  expect_equal(out$all.err$bonferroni, cbind(t0 - bonferroni[, 1L], bonferroni[, 2L] - t0))
  expect_equal(out$all.err$simultaneous, cbind(t0 - simultaneous[, 1L], simultaneous[, 2L] - t0))
})
