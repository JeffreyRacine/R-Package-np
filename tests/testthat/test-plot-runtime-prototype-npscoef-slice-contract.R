test_that("npscoef fixed no-error plot-data prototype matches current route", {
  proto <- getFromNamespace(".np_plot_proto_npscoef_fixed_none_data", "np")
  withr::local_options(np.messages = FALSE)
  set.seed(2601)

  n <- 70L
  x <- data.frame(x = runif(n, -1, 1))
  z <- data.frame(z = runif(n, -0.5, 1.5))
  y <- (1 + 0.5 * x$x) * sin(z$z) + rnorm(n, sd = 0.15)
  cases <- list(
    lc = list(regtype = "lc"),
    ll = list(regtype = "ll"),
    lp1 = list(regtype = "lp", degree = 1L),
    lp2 = list(regtype = "lp", degree = 2L)
  )

  for (case.name in names(cases)) {
    bw.args <- c(
      list(
        xdat = x,
        zdat = z,
        ydat = y,
        bws = 0.45,
        bandwidth.compute = FALSE,
        bwtype = "fixed"
      ),
      cases[[case.name]]
    )
    bw <- do.call(npscoefbw, bw.args)

    old <- suppressWarnings(plot(
      bw,
      xdat = x,
      zdat = z,
      ydat = y,
      output = "data",
      errors = "none",
      neval = 7L,
      perspective = TRUE
    ))
    candidate <- proto(
      bw,
      xdat = x,
      zdat = z,
      ydat = y,
      neval = 7L
    )
    stages <- proto(
      bw,
      xdat = x,
      zdat = z,
      ydat = y,
      neval = 7L,
      return.stages = TRUE
    )

    expect_named(candidate, names(old), info = case.name)
    expect_s3_class(candidate$r1, "smoothcoefficient")
    expect_named(candidate$r1, names(old$r1), info = case.name)
    expect_equal(candidate$r1$eval, old$r1$eval, info = case.name)
    expect_equal(candidate$r1$mean, old$r1$mean, info = case.name)
    expect_equal(candidate$r1$merr, old$r1$merr, info = case.name)
    expect_equal(candidate$r1$bias, old$r1$bias, info = case.name)
    expect_named(stages, c("state", "target_grid", "evaluator", "intervals", "bootstrap", "plot_data"),
                 info = case.name)
    expect_equal(stages$plot_data, candidate, info = case.name)
    expect_null(stages$intervals, info = case.name)
    expect_null(stages$bootstrap, info = case.name)
    expect_equal(stages$evaluator$mean, old$r1$mean, info = case.name)
  }
})

test_that("npscoef plot prototype fails early outside its explicit fixed surface slice", {
  proto <- getFromNamespace(".np_plot_proto_npscoef_fixed_none_data", "np")
  withr::local_options(np.messages = FALSE)
  set.seed(2602)

  n <- 55L
  x <- data.frame(x = runif(n))
  z <- data.frame(z = factor(rbinom(n, 1L, 0.5)))
  y <- x$x + rnorm(n, sd = 0.1)
  bw <- npscoefbw(
    xdat = x,
    zdat = z,
    ydat = y,
    bws = 0.45,
    bandwidth.compute = FALSE,
    regtype = "lc",
    bwtype = "fixed"
  )

  expect_error(
    proto(bw, xdat = x, zdat = z, ydat = y, neval = 7L),
    "continuous/ordered surface variables",
    fixed = TRUE
  )

  expect_error(
    proto(bw, neval = 7L),
    "explicit xdat, ydat, and zdat",
    fixed = TRUE
  )
})

test_that("npscoef fixed asymptotic plot-data prototype matches current route", {
  proto <- getFromNamespace(".np_plot_proto_npscoef_fixed_asymptotic_data", "np")
  withr::local_options(np.messages = FALSE)
  set.seed(2603)

  n <- 65L
  x <- data.frame(x = runif(n, -1, 1))
  z <- data.frame(z = runif(n, -1, 1))
  y <- (1 + x$x) * cos(z$z) + rnorm(n, sd = 0.2)
  bw <- npscoefbw(
    xdat = x,
    zdat = z,
    ydat = y,
    bws = 0.55,
    bandwidth.compute = FALSE,
    regtype = "ll",
    bwtype = "fixed"
  )

  for (band in c("pmzsd", "pointwise", "bonferroni", "simultaneous", "all")) {
    old <- suppressWarnings(plot(
      bw,
      xdat = x,
      zdat = z,
      ydat = y,
      output = "data",
      errors = "asymptotic",
      band = band,
      neval = 6L,
      perspective = TRUE
    ))
    candidate <- proto(
      bw,
      xdat = x,
      zdat = z,
      ydat = y,
      neval = 6L,
      band = band
    )
    stages <- proto(
      bw,
      xdat = x,
      zdat = z,
      ydat = y,
      neval = 6L,
      band = band,
      return.stages = TRUE
    )

    expect_equal(candidate$r1$eval, old$r1$eval, info = band)
    expect_equal(candidate$r1$mean, old$r1$mean, info = band)
    expect_equal(candidate$r1$merr, old$r1$merr, info = band)
    expect_equal(candidate$r1$bias, old$r1$bias, info = band)
    expect_identical(stages$intervals$method, "asymptotic")
    expect_identical(stages$intervals$type, band)
    expect_null(stages$bootstrap, info = band)
    expect_equal(stages$plot_data, candidate, info = band)
  }
})

test_that("npscoef fixed bootstrap plot-data prototype matches current route", {
  proto <- getFromNamespace(".np_plot_proto_npscoef_fixed_bootstrap_data", "np")
  withr::local_options(np.messages = FALSE)
  set.seed(2604)

  n <- 50L
  x <- data.frame(x = runif(n, -1, 1))
  z <- data.frame(z = runif(n, -1, 1))
  y <- (1 - 0.5 * x$x) * sin(z$z) + rnorm(n, sd = 0.25)
  bw <- npscoefbw(
    xdat = x,
    zdat = z,
    ydat = y,
    bws = 0.55,
    bandwidth.compute = FALSE,
    regtype = "ll",
    bwtype = "fixed"
  )
  cases <- expand.grid(
    method = c("wild", "inid"),
    center = c("estimate", "bias-corrected"),
    stringsAsFactors = FALSE
  )

  for (ii in seq_len(nrow(cases))) {
    method <- cases$method[ii]
    center <- cases$center[ii]
    boot.seed <- 6200L + ii
    set.seed(boot.seed)
    old <- suppressWarnings(plot(
      bw,
      xdat = x,
      zdat = z,
      ydat = y,
      output = "data",
      errors = "bootstrap",
      bootstrap = method,
      B = 11L,
      center = center,
      band = "pointwise",
      neval = 5L,
      perspective = TRUE,
      random.seed = boot.seed
    ))
    set.seed(boot.seed)
    candidate <- suppressWarnings(proto(
      bw,
      xdat = x,
      zdat = z,
      ydat = y,
      neval = 5L,
      bootstrap = method,
      B = 11L,
      center = center,
      band = "pointwise"
    ))
    set.seed(boot.seed)
    stages <- suppressWarnings(proto(
      bw,
      xdat = x,
      zdat = z,
      ydat = y,
      neval = 5L,
      bootstrap = method,
      B = 11L,
      center = center,
      band = "pointwise",
      return.stages = TRUE
    ))

    expect_equal(candidate$r1$eval, old$r1$eval, info = paste(method, center))
    expect_equal(candidate$r1$mean, old$r1$mean, info = paste(method, center))
    expect_equal(candidate$r1$merr, old$r1$merr, info = paste(method, center))
    expect_equal(candidate$r1$bias, old$r1$bias, info = paste(method, center))
    expect_identical(stages$bootstrap$method, method)
    expect_identical(stages$bootstrap$center, center)
    expect_equal(stages$plot_data, candidate, info = paste(method, center))
  }
})
