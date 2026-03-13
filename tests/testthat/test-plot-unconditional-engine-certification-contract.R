with_plot_contract_device <- function(expr) {
  pdf(file = tempfile(fileext = ".pdf"))
  on.exit(dev.off(), add = TRUE)
  suppressWarnings(force(expr))
}

expect_plot_modes_match_unconditional <- function(bw,
                                                  xdat,
                                                  method_label,
                                                  yfield,
                                                  plot.errors.method,
                                                  plot.errors.type) {
  data.out <- suppressWarnings(
    plot(
      bw,
      xdat = xdat,
      plot.behavior = "data",
      perspective = FALSE,
      plot.errors.method = plot.errors.method,
      plot.errors.type = plot.errors.type,
      plot.errors.boot.num = 5,
      view = "fixed"
    )
  )

  plot.out <- with_plot_contract_device(
    plot(
      bw,
      xdat = xdat,
      plot.behavior = "plot-data",
      perspective = FALSE,
      plot.errors.method = plot.errors.method,
      plot.errors.type = plot.errors.type,
      plot.errors.boot.num = 5,
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
  expect_true(is.list(state$oldpar))
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
    plot.errors.method = "bootstrap",
    plot.errors.type = "all"
  )
  expect_plot_modes_match_unconditional(
    bw = dist.bw,
    xdat = xdat,
    method_label = "npudistbw bootstrap",
    yfield = "dist",
    plot.errors.method = "bootstrap",
    plot.errors.type = "all"
  )
  expect_plot_modes_match_unconditional(
    bw = dens.bw,
    xdat = xdat,
    method_label = "npudensbw asymptotic",
    yfield = "dens",
    plot.errors.method = "asymptotic",
    plot.errors.type = "pointwise"
  )
  expect_plot_modes_match_unconditional(
    bw = dist.bw,
    xdat = xdat,
    method_label = "npudistbw asymptotic",
    yfield = "dist",
    plot.errors.method = "asymptotic",
    plot.errors.type = "pointwise"
  )
})
