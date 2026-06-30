test_that("npreg fixed no-error plot-data prototype matches current route", {
  proto <- getFromNamespace(".np_plot_proto_npreg_fixed_none_data", "np")
  withr::local_options(np.messages = FALSE)
  set.seed(2201)

  n <- 70L
  x <- data.frame(x1 = runif(n, -1, 1), x2 = runif(n, -0.5, 1.5))
  y <- sin(2 * pi * x$x1) + 0.5 * x$x2 + rnorm(n, sd = 0.15)
  cases <- list(
    lc = list(regtype = "lc"),
    ll = list(regtype = "ll"),
    lp1 = list(regtype = "lp", degree = c(1L, 1L)),
    lp2 = list(regtype = "lp", degree = c(2L, 2L))
  )

  for (case.name in names(cases)) {
    bw.args <- c(
      list(
        xdat = x,
        ydat = y,
        bws = c(0.45, 0.55),
        bandwidth.compute = FALSE,
        bwtype = "fixed"
      ),
      cases[[case.name]]
    )
    bw <- do.call(npregbw, bw.args)

    old <- suppressWarnings(plot(
      bw,
      xdat = x,
      ydat = y,
      output = "data",
      errors = "none",
      view = "fixed",
      neval = 7L,
      perspective = TRUE
    ))
    candidate <- proto(
      bw,
      xdat = x,
      ydat = y,
      neval = 7L
    )
    stages <- proto(
      bw,
      xdat = x,
      ydat = y,
      neval = 7L,
      return.stages = TRUE
    )

    expect_named(candidate, names(old), info = case.name)
    expect_s3_class(candidate$r1, "npregression")
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

test_that("npreg plot prototype fails early outside its fixed 2D slice", {
  proto <- getFromNamespace(".np_plot_proto_npreg_fixed_none_data", "np")
  withr::local_options(np.messages = FALSE)
  set.seed(2202)

  n <- 50L
  x <- data.frame(x1 = runif(n), x2 = factor(rbinom(n, 1L, 0.5)))
  y <- x$x1 + rnorm(n, sd = 0.1)
  bw <- npregbw(
    xdat = x,
    ydat = y,
    bws = c(0.4, 0.2),
    bandwidth.compute = FALSE,
    regtype = "lc",
    bwtype = "fixed"
  )

  expect_error(
    proto(bw, xdat = x, ydat = y, neval = 7L),
    "continuous/ordered surface variables",
    fixed = TRUE
  )

  expect_error(
    proto(bw, neval = 7L),
    "explicit xdat and ydat",
    fixed = TRUE
  )
})

test_that("npreg fixed asymptotic plot-data prototype matches current route", {
  proto <- getFromNamespace(".np_plot_proto_npreg_fixed_asymptotic_data", "np")
  withr::local_options(np.messages = FALSE)
  set.seed(2203)

  n <- 65L
  x <- data.frame(x1 = runif(n, -1, 1), x2 = runif(n, -1, 1))
  y <- x$x1^2 - x$x2 + rnorm(n, sd = 0.2)
  bw <- npregbw(
    xdat = x,
    ydat = y,
    bws = c(0.5, 0.6),
    bandwidth.compute = FALSE,
    regtype = "lp",
    degree = c(1L, 2L),
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
      view = "fixed",
      neval = 6L,
      perspective = TRUE
    ))
    candidate <- proto(
      bw,
      xdat = x,
      ydat = y,
      neval = 6L,
      band = band
    )
    stages <- proto(
      bw,
      xdat = x,
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

test_that("npreg fixed bootstrap plot-data prototype matches current route", {
  proto <- getFromNamespace(".np_plot_proto_npreg_fixed_bootstrap_data", "np")
  withr::local_options(np.messages = FALSE)
  set.seed(2204)

  n <- 55L
  x <- data.frame(x1 = runif(n, -1, 1), x2 = runif(n, -1, 1))
  y <- cos(2 * x$x1) + x$x2 + rnorm(n, sd = 0.25)
  bw <- npregbw(
    xdat = x,
    ydat = y,
    bws = c(0.55, 0.65),
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
    boot.seed <- 3000L + ii
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
      view = "fixed",
      neval = 5L,
      perspective = TRUE,
      random.seed = boot.seed
    ))
    set.seed(boot.seed)
    candidate <- suppressWarnings(proto(
      bw,
      xdat = x,
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
      ydat = y,
      neval = 5L,
      bootstrap = method,
      B = 11L,
      center = center,
      band = "pointwise",
      return.stages = TRUE
    ))

    label <- paste(method, center)
    expect_equal(candidate$r1$eval, old$r1$eval, info = label)
    expect_equal(candidate$r1$mean, old$r1$mean, info = label)
    expect_equal(candidate$r1$merr, old$r1$merr, info = label)
    expect_equal(candidate$r1$bias, old$r1$bias, info = label)
    expect_identical(stages$intervals$method, "bootstrap")
    expect_identical(stages$intervals$type, "pointwise")
    expect_identical(stages$bootstrap$method, method)
    expect_identical(stages$bootstrap$B, 11L)
    expect_equal(stages$plot_data, candidate, info = label)
  }
})

test_that("npreg staged plot data can be rendered without estimator re-entry", {
  proto <- getFromNamespace(".np_plot_proto_npreg_fixed_none_data", "np")
  render <- getFromNamespace(".np_plot_proto_rectangular_surface_base_render", "np")
  withr::local_options(np.messages = FALSE)
  set.seed(2205)

  n <- 55L
  x <- data.frame(x1 = runif(n, -1, 1), x2 = runif(n, -1, 1))
  y <- x$x1 + x$x2^2 + rnorm(n, sd = 0.1)
  bw <- npregbw(
    xdat = x,
    ydat = y,
    bws = c(0.55, 0.55),
    bandwidth.compute = FALSE,
    regtype = "lp",
    degree = c(1L, 1L),
    bwtype = "fixed"
  )
  plot.data <- proto(bw, xdat = x, ydat = y, neval = 5L)
  pdf.file <- tempfile(fileext = ".pdf")
  grDevices::pdf(pdf.file)
  on.exit(if (grDevices::dev.cur() > 1L) grDevices::dev.off(), add = TRUE)

  out.image <- render(plot.data, perspective = FALSE)
  out.persp <- render(plot.data, perspective = TRUE)

  expect_identical(out.image, plot.data)
  expect_identical(out.persp, plot.data)
  grDevices::dev.off()
  expect_true(file.exists(pdf.file))
  expect_gt(file.info(pdf.file)$size, 0)
})
