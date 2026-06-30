test_that("npindex fixed no-error plot-data prototype matches current route", {
  proto <- getFromNamespace(".np_plot_proto_npindex_fixed_none_data", "np")
  withr::local_options(np.messages = FALSE)
  set.seed(2501)

  n <- 70L
  x <- data.frame(x1 = runif(n, -1, 1), x2 = runif(n, -0.5, 1.5))
  y <- sin(x$x1 + x$x2) + rnorm(n, sd = 0.15)
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
        ydat = y,
        bws = c(1, 1, 0.45),
        bandwidth.compute = FALSE,
        bwtype = "fixed"
      ),
      cases[[case.name]]
    )
    bw <- do.call(npindexbw, bw.args)

    old <- suppressWarnings(plot(
      bw,
      xdat = x,
      ydat = y,
      output = "data",
      errors = "none",
      neval = 17L,
      gradients = FALSE
    ))
    candidate <- proto(
      bw,
      xdat = x,
      ydat = y,
      neval = 17L
    )
    stages <- proto(
      bw,
      xdat = x,
      ydat = y,
      neval = 17L,
      return.stages = TRUE
    )

    expect_named(candidate, names(old), info = case.name)
    expect_s3_class(candidate$si1, "singleindex")
    expect_named(candidate$si1, names(old$si1), info = case.name)
    expect_equal(candidate$si1$index, old$si1$index, info = case.name)
    expect_equal(candidate$si1$mean, old$si1$mean, info = case.name)
    expect_equal(candidate$si1$merr, old$si1$merr, info = case.name)
    expect_equal(candidate$si1$bias, old$si1$bias, info = case.name)
    expect_false(candidate$si1$gradients, info = case.name)
    expect_named(stages, c("state", "target_grid", "evaluator", "intervals", "bootstrap", "plot_data"),
                 info = case.name)
    expect_equal(stages$plot_data, candidate, info = case.name)
    expect_null(stages$intervals, info = case.name)
    expect_null(stages$bootstrap, info = case.name)
    expect_equal(stages$evaluator$mean, old$si1$mean, info = case.name)
  }
})

test_that("npindex plot prototype fails early outside its fixed non-gradient slice", {
  proto <- getFromNamespace(".np_plot_proto_npindex_fixed_none_data", "np")
  withr::local_options(np.messages = FALSE)
  set.seed(2502)

  n <- 55L
  x <- data.frame(x1 = runif(n), x2 = factor(rbinom(n, 1L, 0.5)))
  y <- x$x1 + rnorm(n, sd = 0.1)
  bw <- npindexbw(
    xdat = x,
    ydat = y,
    bws = c(1, 1, 0.45),
    bandwidth.compute = FALSE,
    regtype = "lc",
    bwtype = "fixed"
  )

  expect_error(
    proto(bw, xdat = x, ydat = y, neval = 17L),
    "continuous/ordered index variables",
    fixed = TRUE
  )

  expect_error(
    proto(bw, neval = 17L),
    "explicit xdat and ydat",
    fixed = TRUE
  )
})

test_that("npindex fixed asymptotic plot-data prototype matches current route", {
  proto <- getFromNamespace(".np_plot_proto_npindex_fixed_asymptotic_data", "np")
  withr::local_options(np.messages = FALSE)
  set.seed(2503)

  n <- 65L
  x <- data.frame(x1 = runif(n, -1, 1), x2 = runif(n, -1, 1))
  y <- x$x1^2 - x$x2 + rnorm(n, sd = 0.2)
  bw <- npindexbw(
    xdat = x,
    ydat = y,
    bws = c(1, 1, 0.55),
    bandwidth.compute = FALSE,
    regtype = "lp",
    degree = 2L,
    bwtype = "fixed"
  )

  for (band in c("pmzsd", "pointwise", "bonferroni", "simultaneous", "all")) {
    old <- suppressWarnings(plot(
      bw,
      xdat = x,
      ydat = y,
      output = "data",
      errors = "asymptotic",
      band = band,
      neval = 15L,
      gradients = FALSE
    ))
    candidate <- proto(
      bw,
      xdat = x,
      ydat = y,
      neval = 15L,
      band = band
    )
    stages <- proto(
      bw,
      xdat = x,
      ydat = y,
      neval = 15L,
      band = band,
      return.stages = TRUE
    )

    expect_equal(candidate$si1$index, old$si1$index, info = band)
    expect_equal(candidate$si1$mean, old$si1$mean, info = band)
    expect_equal(candidate$si1$merr, old$si1$merr, info = band)
    expect_equal(candidate$si1$bias, old$si1$bias, info = band)
    expect_identical(stages$intervals$method, "asymptotic")
    expect_identical(stages$intervals$type, band)
    expect_null(stages$bootstrap, info = band)
    expect_equal(stages$plot_data, candidate, info = band)
  }
})

test_that("npindex fixed bootstrap plot-data prototype matches current route", {
  proto <- getFromNamespace(".np_plot_proto_npindex_fixed_bootstrap_data", "np")
  withr::local_options(np.messages = FALSE)
  set.seed(2504)

  n <- 55L
  x <- data.frame(x1 = runif(n, -1, 1), x2 = runif(n, -1, 1))
  y <- cos(x$x1 + x$x2) + rnorm(n, sd = 0.25)
  bw <- npindexbw(
    xdat = x,
    ydat = y,
    bws = c(1, 1, 0.55),
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
    boot.seed <- 5200L + ii
    set.seed(boot.seed)
    old <- suppressWarnings(plot(
      bw,
      xdat = x,
      ydat = y,
      output = "data",
      errors = "bootstrap",
      bootstrap = method,
      B = 11L,
      center = center,
      band = "pointwise",
      neval = 13L,
      gradients = FALSE,
      random.seed = boot.seed
    ))
    set.seed(boot.seed)
    candidate <- suppressWarnings(proto(
      bw,
      xdat = x,
      ydat = y,
      neval = 13L,
      bootstrap = method,
      B = 11L,
      center = center,
      band = "pointwise"
    ))
    set.seed(boot.seed)
    stages <- suppressWarnings(proto(
      bw,
      xdat = x,
      ydat = y,
      neval = 13L,
      bootstrap = method,
      B = 11L,
      center = center,
      band = "pointwise",
      return.stages = TRUE
    ))

    expect_equal(candidate$si1$index, old$si1$index, info = paste(method, center))
    expect_equal(candidate$si1$mean, old$si1$mean, info = paste(method, center))
    expect_equal(candidate$si1$merr, old$si1$merr, info = paste(method, center))
    expect_equal(candidate$si1$bias, old$si1$bias, info = paste(method, center))
    expect_identical(stages$bootstrap$method, method)
    expect_identical(stages$bootstrap$center, center)
    expect_equal(stages$plot_data, candidate, info = paste(method, center))
  }
})
