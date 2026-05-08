with_plot_contract_device <- function(expr) {
  pdf(file = tempfile(fileext = ".pdf"))
  on.exit(dev.off(), add = TRUE)
  suppressWarnings(force(expr))
}

expect_plot_modes_match_unconditional <- function(bw,
                                                  xdat,
                                                  method_label,
                                                  yfield,
                                                  errors,
                                                  band) {
  data.out <- suppressWarnings(
    plot(
      bw,
      xdat = xdat,
      output = "data",
      perspective = FALSE,
      errors = errors,
      band = band,
      B = 5,
      view = "fixed"
    )
  )

  plot.out <- with_plot_contract_device(
    plot(
      bw,
      xdat = xdat,
      output = "plot-data",
      perspective = FALSE,
      errors = errors,
      band = band,
      B = 5,
      view = "fixed"
    )
  )

  expect_true(is.list(plot.out), info = method_label)
  expect_named(plot.out, names(data.out), info = method_label)

  for (nm in names(data.out)) {
    expect_true(
      inherits(plot.out[[nm]], class(data.out[[nm]])[1L]),
      info = paste(method_label, nm, "class")
    )
    expect_equal(plot.out[[nm]]$eval, data.out[[nm]]$eval, info = paste(method_label, nm, "eval"))
    expect_equal(plot.out[[nm]][[yfield]], data.out[[nm]][[yfield]], info = paste(method_label, nm, yfield))
    expect_equal(plot.out[[nm]]$derr, data.out[[nm]]$derr, info = paste(method_label, nm, "derr"))
    expect_equal(plot.out[[nm]]$bias, data.out[[nm]]$bias, info = paste(method_label, nm, "bias"))
  }
}

test_that("plot engine startup helper respects plot.par.mfrow option override", {
  skip_if_not_installed("np")

  begin <- getFromNamespace(".np_plot_engine_begin", "np")
  restore <- getFromNamespace(".np_plot_restore_par", "np")

  pdf(file = tempfile(fileext = ".pdf"))
  on.exit(dev.off(), add = TRUE)

  old_opts <- options(plot.par.mfrow = FALSE)
  on.exit(options(old_opts), add = TRUE)

  state <- begin(plot.par.mfrow = TRUE)
  on.exit(restore(state$oldpar), add = TRUE)

  expect_identical(state$plot.par.mfrow, FALSE)
  expect_true(is.numeric(state$oldpar))
  expect_length(state$oldpar, 1L)
})

test_that("unconditional interval payload helper preserves bootstrap and asymptotic semantics", {
  skip_if_not_installed("np")

  payload <- getFromNamespace(".np_plot_interval_payload", "np")

  boot.raw <- list(
    boot.err = cbind(c(0.2, 0.3), c(0.4, 0.5), c(1.2, 1.8)),
    boot.all.err = list(pointwise = cbind(c(0.2, 0.3), c(0.4, 0.5))),
    bxp = list(stats = matrix(1, nrow = 5, ncol = 1))
  )

  boot.out <- payload(
    estimate = c(1, 2),
    se = c(0.1, 0.2),
    errors = "bootstrap",
    alpha = 0.05,
    band = "all",
    center = "bias-corrected",
    bootstrap_raw = boot.raw
  )

  expect_equal(boot.out$err, boot.raw$boot.err)
  expect_identical(boot.out$all.err, boot.raw$boot.all.err)
  expect_identical(boot.out$center, boot.raw$boot.err[,3])
  expect_identical(boot.out$bxp, boot.raw$bxp)

  asym.out <- payload(
    estimate = c(1, 2),
    se = c(0.1, 0.2),
    errors = "asymptotic",
    alpha = 0.05,
    band = "pointwise",
    center = "estimate",
    bootstrap_raw = NULL
  )

  expect_identical(asym.out$center, c(1, 2))
  expect_equal(dim(asym.out$err), c(2L, 3L))
  expect_true(all(is.na(asym.out$err[,3])))
  expect_identical(asym.out$bxp, list())
})

test_that("bandwidth and dbandwidth engines preserve data vs plot-data payloads", {
  skip_if_not_installed("np")

  set.seed(20260313)
  xdat <- data.frame(x = rnorm(18), z = runif(18))

  dens.bw <- npudensbw(
    dat = xdat,
    bws = c(0.35, 0.35),
    bandwidth.compute = FALSE
  )
  dist.bw <- npudistbw(
    dat = xdat,
    bws = c(0.35, 0.35),
    bandwidth.compute = FALSE
  )

  expect_plot_modes_match_unconditional(
    bw = dens.bw,
    xdat = xdat,
    method_label = "npudensbw bootstrap",
    yfield = "dens",
    errors = "bootstrap",
    band = "all"
  )
  expect_plot_modes_match_unconditional(
    bw = dist.bw,
    xdat = xdat,
    method_label = "npudistbw bootstrap",
    yfield = "dist",
    errors = "bootstrap",
    band = "all"
  )
  expect_plot_modes_match_unconditional(
    bw = dens.bw,
    xdat = xdat,
    method_label = "npudensbw asymptotic",
    yfield = "dens",
    errors = "asymptotic",
    band = "pointwise"
  )
  expect_plot_modes_match_unconditional(
    bw = dist.bw,
    xdat = xdat,
    method_label = "npudistbw asymptotic",
    yfield = "dist",
    errors = "asymptotic",
    band = "pointwise"
  )
})
