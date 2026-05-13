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
    output = "data",
    errors = "none",
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
      output = "data",
      errors = "asymptotic",
      band = band,
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

test_that("npqreg fixed quantile-level bootstrap slice prototype matches current fitted-object route", {
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

  cases <- expand.grid(
    method = c("inid", "fixed", "geom"),
    center = c("estimate", "bias-corrected"),
    stringsAsFactors = FALSE
  )
  for (ii in seq_len(nrow(cases))) {
    method <- cases$method[ii]
    center <- cases$center[ii]
    boot.seed <- 8200L + ii
    plot.extra <- if (identical(method, "inid")) list() else list(boot_control = np_boot_control(blocklen = 3L))
    proto.extra <- if (identical(method, "inid")) list() else list(plot.errors.boot.blocklen = 3L)
    set.seed(boot.seed)
    old <- suppressWarnings(do.call(plot, c(
      list(
        x = fit,
        xdat = x,
        ydat = y,
        output = "data",
        errors = "bootstrap",
        bootstrap = method,
        B = 11L,
        center = center,
        band = "pmzsd",
        neval = 5L,
        perspective = FALSE,
        random.seed = boot.seed
      ),
      plot.extra
    )))
    set.seed(boot.seed)
    candidate <- suppressWarnings(do.call(proto, c(
      list(
        object = fit,
        xdat = x,
        ydat = y,
        neval = 5L,
        plot.errors.boot.method = method,
        plot.errors.boot.num = 11L,
        plot.errors.center = center,
        plot.errors.type = "pmzsd"
      ),
      proto.extra
    )))
    set.seed(boot.seed)
    stages <- suppressWarnings(do.call(proto, c(
      list(
        object = fit,
        xdat = x,
        ydat = y,
        neval = 5L,
        plot.errors.boot.method = method,
        plot.errors.boot.num = 11L,
        plot.errors.center = center,
        plot.errors.type = "pmzsd",
        return.stages = TRUE
      ),
      proto.extra
    )))

    for (nm in names(old)) {
      expect_equal(candidate[[nm]]$xeval, old[[nm]]$xeval, info = paste(nm, method, center))
      expect_equal(candidate[[nm]]$quantile, old[[nm]]$quantile, info = paste(nm, method, center))
      expect_equal(candidate[[nm]]$quanterr, old[[nm]]$quanterr, info = paste(nm, method, center))
      expect_equal(candidate[[nm]]$bias, old[[nm]]$bias, info = paste(nm, method, center))
      expect_equal(candidate[[nm]]$bxp, old[[nm]]$bxp, info = paste(nm, method, center))
    }
    expect_true(all(vapply(stages$bootstrap, function(x) is.list(x) && !is.null(x$boot.err), logical(1))),
                info = paste(method, center))
  }
})

test_that("npqreg fixed gradient slice prototype matches current fitted-object route", {
  proto.none <- getFromNamespace(".np_plot_proto_npqreg_fixed_none_data", "np")
  proto.boot <- getFromNamespace(".np_plot_proto_npqreg_fixed_bootstrap_data", "np")
  proto.asym <- getFromNamespace(".np_plot_proto_npqreg_fixed_asymptotic_data", "np")
  withr::local_options(np.messages = FALSE)
  set.seed(2805)

  n <- 45L
  x <- data.frame(x1 = runif(n, -1, 1), x2 = runif(n, -1, 1))
  y <- data.frame(y = x$x1^2 - x$x2 + rnorm(n, sd = 0.15))
  bw <- npcdistbw(
    xdat = x,
    ydat = y,
    bws = c(0.50, 0.55, 0.60),
    bandwidth.compute = FALSE,
    bwtype = "fixed"
  )
  fit <- npqreg(
    bws = bw,
    txdat = x,
    tydat = y,
    tau = 0.45
  )

  old.none <- suppressWarnings(plot(
    fit,
    xdat = x,
    ydat = y,
    output = "data",
    errors = "none",
    gradients = TRUE,
    neval = 5L,
    perspective = FALSE
  ))
  candidate.none <- proto.none(
    fit,
    xdat = x,
    ydat = y,
    neval = 5L,
    gradients = TRUE
  )
  stages.none <- proto.none(
    fit,
    xdat = x,
    ydat = y,
    neval = 5L,
    gradients = TRUE,
    return.stages = TRUE
  )

  expect_named(candidate.none, names(old.none))
  expect_true(isTRUE(stages.none$state$gradients))
  for (nm in names(old.none)) {
    expect_equal(candidate.none[[nm]]$xeval, old.none[[nm]]$xeval, info = nm)
    expect_equal(candidate.none[[nm]]$quantile, old.none[[nm]]$quantile, info = nm)
    expect_equal(candidate.none[[nm]]$quantgrad, old.none[[nm]]$quantgrad, info = nm)
    expect_equal(candidate.none[[nm]]$quanterr, old.none[[nm]]$quanterr, info = nm)
  }

  old.asym <- suppressWarnings(plot(
    fit,
    xdat = x,
    ydat = y,
    output = "data",
    errors = "asymptotic",
    gradients = TRUE,
    neval = 5L,
    perspective = FALSE
  ))
  candidate.asym <- proto.asym(
    fit,
    xdat = x,
    ydat = y,
    neval = 5L,
    gradients = TRUE
  )
  expect_named(candidate.asym, names(old.asym))
  for (nm in names(old.asym)) {
    expect_equal(candidate.asym[[nm]]$xeval, old.asym[[nm]]$xeval, info = nm)
    expect_equal(candidate.asym[[nm]]$quantgrad, old.asym[[nm]]$quantgrad, info = nm)
    expect_equal(candidate.asym[[nm]]$quantgerr, old.asym[[nm]]$quantgerr, info = nm)
  }

  for (method in c("inid", "fixed")) {
    boot.seed <- if (identical(method, "inid")) 8301L else 902L
    plot.extra <- if (identical(method, "fixed")) list(boot_control = np_boot_control(blocklen = 2L)) else list()
    proto.extra <- if (identical(method, "fixed")) list(plot.errors.boot.blocklen = 2L) else list()
    boot.num <- if (identical(method, "fixed")) 5L else 7L
    set.seed(boot.seed)
    old.boot <- suppressWarnings(do.call(plot, c(
      list(
        x = fit,
        xdat = x,
        ydat = y,
        output = "data",
        errors = "bootstrap",
        bootstrap = method,
        B = boot.num,
        center = "estimate",
        band = "pmzsd",
        gradients = TRUE,
        neval = 5L,
        perspective = FALSE,
        random.seed = boot.seed
      ),
      plot.extra
    )))
    set.seed(boot.seed)
    candidate.boot <- suppressWarnings(do.call(proto.boot, c(
      list(
        object = fit,
        xdat = x,
        ydat = y,
        neval = 5L,
        plot.errors.boot.method = method,
        plot.errors.boot.num = boot.num,
        plot.errors.center = "estimate",
        plot.errors.type = "pmzsd",
        gradients = TRUE
      ),
      proto.extra
    )))
    set.seed(boot.seed)
    stages.boot <- suppressWarnings(do.call(proto.boot, c(
      list(
        object = fit,
        xdat = x,
        ydat = y,
        neval = 5L,
        plot.errors.boot.method = method,
        plot.errors.boot.num = boot.num,
        plot.errors.center = "estimate",
        plot.errors.type = "pmzsd",
        gradients = TRUE,
        return.stages = TRUE
      ),
      proto.extra
    )))

    for (nm in names(old.boot)) {
      expect_equal(candidate.boot[[nm]]$xeval, old.boot[[nm]]$xeval, info = paste(nm, method))
      expect_equal(candidate.boot[[nm]]$quantile, old.boot[[nm]]$quantile, info = paste(nm, method))
      expect_equal(candidate.boot[[nm]]$quantgrad, old.boot[[nm]]$quantgrad, info = paste(nm, method))
      expect_equal(candidate.boot[[nm]]$gc1err, old.boot[[nm]]$gc1err, info = paste(nm, method))
      expect_equal(candidate.boot[[nm]]$gc1bias, old.boot[[nm]]$gc1bias, info = paste(nm, method))
      expect_equal(candidate.boot[[nm]]$gc2err, old.boot[[nm]]$gc2err, info = paste(nm, method))
      expect_equal(candidate.boot[[nm]]$gc2bias, old.boot[[nm]]$gc2bias, info = paste(nm, method))
      expect_equal(candidate.boot[[nm]]$bxp, old.boot[[nm]]$bxp, info = paste(nm, method))
    }
    expect_true(all(vapply(stages.boot$bootstrap, is.list, logical(1))))
  }
})
