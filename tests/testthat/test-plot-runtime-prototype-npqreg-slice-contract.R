test_that("npqreg fixed no-error slice prototype matches current fitted-object route", {
  proto <- getFromNamespace(".np_plot_proto_npqreg_fixed_none_data", "np")
  withr::local_options(np.messages = FALSE)
  set.seed(2801)

  n <- 70L
  x <- data.frame(x1 = runif(n, -1, 1), x2 = runif(n, -0.5, 1.5))
  y <- data.frame(y = x$x1 - x$x2 + rnorm(n, sd = 0.15))
  bw <- npcdistbw(
    xdat = x,
    ydat = y,
    bws = c(0.45, 0.50, 0.55),
    bandwidth.compute = FALSE,
    bwtype = "fixed"
  )
  fit <- npqreg(
    bws = bw,
    txdat = x,
    tydat = y,
    tau = 0.4
  )

  old <- suppressWarnings(plot(
    fit,
    xdat = x,
    ydat = y,
    plot.behavior = "data",
    plot.errors.method = "none",
    neval = 7L,
    perspective = FALSE
  ))
  candidate <- proto(
    fit,
    xdat = x,
    ydat = y,
    neval = 7L
  )
  stages <- proto(
    fit,
    xdat = x,
    ydat = y,
    neval = 7L,
    return.stages = TRUE
  )

  expect_named(candidate, names(old))
  expect_true(all(vapply(candidate, inherits, logical(1), "qregression")))
  for (nm in names(old)) {
    expect_named(candidate[[nm]], names(old[[nm]]), info = nm)
    expect_equal(candidate[[nm]]$xeval, old[[nm]]$xeval, info = nm)
    expect_equal(candidate[[nm]]$tau, old[[nm]]$tau, info = nm)
    expect_equal(candidate[[nm]]$quantile, old[[nm]]$quantile, info = nm)
    expect_equal(candidate[[nm]]$quanterr, old[[nm]]$quanterr, info = nm)
  }
  expect_named(stages, c("state", "target_grid", "evaluator", "intervals", "bootstrap", "plot_data"))
  for (nm in names(candidate)) {
    expect_equal(stages$plot_data[[nm]]$xeval, candidate[[nm]]$xeval, info = nm)
    expect_equal(stages$plot_data[[nm]]$quantile, candidate[[nm]]$quantile, info = nm)
    expect_equal(stages$plot_data[[nm]]$quanterr, candidate[[nm]]$quanterr, info = nm)
  }
  expect_true(all(vapply(stages$intervals, is.null, logical(1))))
  expect_true(all(vapply(stages$bootstrap, is.null, logical(1))))
})

test_that("npqreg plot prototype fails early outside its fitted-object fixed slice", {
  proto <- getFromNamespace(".np_plot_proto_npqreg_fixed_none_data", "np")
  withr::local_options(np.messages = FALSE)
  set.seed(2802)

  n <- 55L
  x <- data.frame(x1 = runif(n), x2 = factor(rbinom(n, 1L, 0.5)))
  y <- data.frame(y = x$x1 + rnorm(n, sd = 0.1))
  bw <- npcdistbw(
    xdat = x,
    ydat = y,
    bws = c(0.45, 0.49, 0.49),
    bandwidth.compute = FALSE,
    bwtype = "fixed"
  )
  fit <- npqreg(
    bws = bw,
    txdat = x,
    tydat = y,
    tau = 0.4
  )

  expect_error(
    proto(bw, xdat = x, ydat = y, neval = 7L),
    "quantile regression object",
    fixed = TRUE
  )

  expect_error(
    proto(fit, xdat = x, ydat = y, neval = 7L),
    "continuous/ordered slice variables",
    fixed = TRUE
  )

  expect_error(
    proto(fit, neval = 7L),
    "explicit xdat and ydat",
    fixed = TRUE
  )
})

test_that("npqreg fixed asymptotic slice prototype matches current fitted-object route", {
  proto <- getFromNamespace(".np_plot_proto_npqreg_fixed_asymptotic_data", "np")
  withr::local_options(np.messages = FALSE)
  set.seed(2803)

  n <- 65L
  x <- data.frame(x1 = runif(n, -1, 1), x2 = runif(n, -1, 1))
  y <- data.frame(y = x$x1^2 - x$x2 + rnorm(n, sd = 0.2))
  bw <- npcdistbw(
    xdat = x,
    ydat = y,
    bws = c(0.55, 0.50, 0.60),
    bandwidth.compute = FALSE,
    bwtype = "fixed"
  )
  fit <- npqreg(
    bws = bw,
    txdat = x,
    tydat = y,
    tau = 0.35
  )

  for (band in c("pmzsd", "pointwise", "bonferroni", "simultaneous", "all")) {
    old <- suppressWarnings(plot(
      fit,
      xdat = x,
      ydat = y,
      plot.behavior = "data",
      plot.errors.method = "asymptotic",
      plot.errors.type = band,
      neval = 6L,
      perspective = FALSE
    ))
    candidate <- proto(
      fit,
      xdat = x,
      ydat = y,
      neval = 6L,
      plot.errors.type = band
    )
    stages <- proto(
      fit,
      xdat = x,
      ydat = y,
      neval = 6L,
      plot.errors.type = band,
      return.stages = TRUE
    )

    for (nm in names(old)) {
      expect_equal(candidate[[nm]]$xeval, old[[nm]]$xeval, info = paste(nm, band))
      expect_equal(candidate[[nm]]$quantile, old[[nm]]$quantile, info = paste(nm, band))
      expect_equal(candidate[[nm]]$quanterr, old[[nm]]$quanterr, info = paste(nm, band))
      expect_equal(candidate[[nm]]$bias, old[[nm]]$bias, info = paste(nm, band))
      expect_equal(candidate[[nm]]$bxp, old[[nm]]$bxp, info = paste(nm, band))
    }
    expect_true(all(vapply(stages$bootstrap, is.null, logical(1))))
    for (nm in names(candidate)) {
      expect_equal(stages$plot_data[[nm]]$xeval, candidate[[nm]]$xeval, info = paste(nm, band))
      expect_equal(stages$plot_data[[nm]]$quantile, candidate[[nm]]$quantile, info = paste(nm, band))
      expect_equal(stages$plot_data[[nm]]$quanterr, candidate[[nm]]$quanterr, info = paste(nm, band))
      expect_equal(stages$plot_data[[nm]]$bias, candidate[[nm]]$bias, info = paste(nm, band))
      expect_equal(stages$plot_data[[nm]]$bxp, candidate[[nm]]$bxp, info = paste(nm, band))
    }
  }
})

test_that("npqreg fixed bootstrap slice prototype preserves current fail-fast boundary", {
  proto <- getFromNamespace(".np_plot_proto_npqreg_fixed_bootstrap_data", "np")
  withr::local_options(np.messages = FALSE)
  set.seed(2804)

  n <- 50L
  x <- data.frame(x1 = runif(n, -1, 1), x2 = runif(n, -1, 1))
  y <- data.frame(y = cos(x$x1) + x$x2 + rnorm(n, sd = 0.25))
  bw <- npcdistbw(
    xdat = x,
    ydat = y,
    bws = c(0.55, 0.60, 0.65),
    bandwidth.compute = FALSE,
    bwtype = "fixed"
  )
  fit <- npqreg(
    bws = bw,
    txdat = x,
    tydat = y,
    tau = 0.45
  )

  for (center in c("estimate", "bias-corrected")) {
    boot.seed <- if (identical(center, "estimate")) 8201L else 8202L
    set.seed(boot.seed)
    expect_error(
      suppressWarnings(plot(
        fit,
        xdat = x,
        ydat = y,
        plot.behavior = "data",
        plot.errors.method = "bootstrap",
        plot.errors.boot.method = "inid",
        plot.errors.boot.num = 11L,
        plot.errors.center = center,
        plot.errors.type = "pointwise",
        neval = 5L,
        perspective = FALSE,
        random.seed = boot.seed
      )),
      "inid conditional helper unavailable",
      fixed = TRUE
    )
    set.seed(boot.seed)
    expect_error(
      suppressWarnings(proto(
        fit,
        xdat = x,
        ydat = y,
        neval = 5L,
        plot.errors.boot.method = "inid",
        plot.errors.boot.num = 11L,
        plot.errors.center = center,
        plot.errors.type = "pointwise"
      )),
      "inid conditional helper unavailable",
      fixed = TRUE
    )
  }
})
