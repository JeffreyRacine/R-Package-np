test_that("npudens and npudist fixed no-error plot-data prototypes match current routes", {
  dens.proto <- getFromNamespace(".np_plot_proto_npudens_fixed_none_data", "np")
  dist.proto <- getFromNamespace(".np_plot_proto_npudist_fixed_none_data", "np")
  withr::local_options(np.messages = FALSE)
  set.seed(2401)

  n <- 75L
  x <- data.frame(x1 = rnorm(n), x2 = runif(n, -1, 1))
  dens.bw <- npudensbw(
    dat = x,
    bws = c(0.55, 0.60),
    bandwidth.compute = FALSE,
    bwtype = "fixed"
  )
  dist.bw <- npudistbw(
    dat = x,
    bws = c(0.55, 0.60),
    bandwidth.compute = FALSE,
    bwtype = "fixed"
  )

  old.dens <- suppressWarnings(plot(
    dens.bw,
    xdat = x,
    output = "data",
    errors = "none",
    view = "fixed",
    neval = 7L,
    perspective = TRUE
  ))
  cand.dens <- dens.proto(
    dens.bw,
    xdat = x,
    neval = 7L
  )
  stages.dens <- dens.proto(
    dens.bw,
    xdat = x,
    neval = 7L,
    return.stages = TRUE
  )

  expect_named(cand.dens, names(old.dens))
  expect_s3_class(cand.dens$d1, "npdensity")
  expect_named(cand.dens$d1, names(old.dens$d1))
  expect_equal(cand.dens$d1$eval, old.dens$d1$eval)
  expect_equal(cand.dens$d1$dens, old.dens$d1$dens)
  expect_equal(cand.dens$d1$derr, old.dens$d1$derr)
  expect_equal(cand.dens$d1$bias, old.dens$d1$bias)
  expect_named(stages.dens, c("state", "target_grid", "evaluator", "intervals", "bootstrap", "plot_data"))
  expect_equal(stages.dens$plot_data, cand.dens)
  expect_null(stages.dens$intervals)
  expect_null(stages.dens$bootstrap)
  expect_equal(stages.dens$evaluator$dens, old.dens$d1$dens)

  old.dist <- suppressWarnings(plot(
    dist.bw,
    xdat = x,
    output = "data",
    errors = "none",
    view = "fixed",
    neval = 7L,
    perspective = TRUE
  ))
  cand.dist <- dist.proto(
    dist.bw,
    xdat = x,
    neval = 7L
  )
  stages.dist <- dist.proto(
    dist.bw,
    xdat = x,
    neval = 7L,
    return.stages = TRUE
  )

  expect_named(cand.dist, names(old.dist))
  expect_s3_class(cand.dist$d1, "npdistribution")
  expect_named(cand.dist$d1, names(old.dist$d1))
  expect_equal(cand.dist$d1$eval, old.dist$d1$eval)
  expect_equal(cand.dist$d1$dist, old.dist$d1$dist)
  expect_equal(cand.dist$d1$derr, old.dist$d1$derr)
  expect_equal(cand.dist$d1$bias, old.dist$d1$bias)
  expect_named(stages.dist, c("state", "target_grid", "evaluator", "intervals", "bootstrap", "plot_data"))
  expect_equal(stages.dist$plot_data, cand.dist)
  expect_null(stages.dist$intervals)
  expect_null(stages.dist$bootstrap)
  expect_equal(stages.dist$evaluator$dist, old.dist$d1$dist)
})

test_that("npudens and npudist fixed asymptotic and bootstrap prototypes match current routes", {
  proto <- getFromNamespace(".np_plot_proto_unconditional_fixed_data", "np")
  withr::local_options(np.messages = FALSE)
  set.seed(2402)

  n <- 65L
  x <- data.frame(x1 = rnorm(n), x2 = runif(n, -1, 1))
  dens.bw <- npudensbw(
    dat = x,
    bws = c(0.55, 0.60),
    bandwidth.compute = FALSE,
    bwtype = "fixed"
  )
  dist.bw <- npudistbw(
    dat = x,
    bws = c(0.55, 0.60),
    bandwidth.compute = FALSE,
    bwtype = "fixed"
  )

  for (case in c("dens", "dist")) {
    bws <- if (identical(case, "dens")) dens.bw else dist.bw
    cdf <- identical(case, "dist")
    value.name <- if (cdf) "dist" else "dens"
    error.name <- if (cdf) "derr" else "derr"

    for (band in c("pmzsd", "pointwise", "bonferroni")) {
      old <- suppressWarnings(plot(
        bws,
        xdat = x,
        output = "data",
        errors = "asymptotic",
        band = band,
        view = "fixed",
        neval = 6L,
        perspective = TRUE
      ))
      candidate <- proto(
        bws,
        xdat = x,
        neval = 6L,
        errors = "asymptotic",
        band = band,
        cdf = cdf
      )
      expect_equal(candidate$d1[[value.name]], old$d1[[value.name]], info = paste(case, band))
      expect_equal(candidate$d1[[error.name]], old$d1[[error.name]], info = paste(case, band))
    }

    boot.seed <- if (cdf) 4202L else 4201L
    set.seed(boot.seed)
    old <- suppressWarnings(plot(
      bws,
      xdat = x,
      output = "data",
      errors = "bootstrap",
      bootstrap = "inid",
      B = 13L,
      center = "bias-corrected",
      band = "pointwise",
      view = "fixed",
      neval = 5L,
      perspective = TRUE,
      random.seed = boot.seed
    ))
    set.seed(boot.seed)
    candidate <- suppressWarnings(proto(
      bws,
      xdat = x,
      neval = 5L,
      errors = "bootstrap",
      bootstrap = "inid",
      B = 13L,
      center = "bias-corrected",
      band = "pointwise",
      cdf = cdf
    ))
    set.seed(boot.seed)
    stages <- suppressWarnings(proto(
      bws,
      xdat = x,
      neval = 5L,
      errors = "bootstrap",
      bootstrap = "inid",
      B = 13L,
      center = "bias-corrected",
      band = "pointwise",
      cdf = cdf,
      return.stages = TRUE
    ))

    expect_equal(candidate$d1[[value.name]], old$d1[[value.name]], info = case)
    expect_equal(candidate$d1[[error.name]], old$d1[[error.name]], info = case)
    expect_equal(candidate$d1$bias, old$d1$bias, info = case)
    expect_identical(stages$intervals$method, "bootstrap")
    expect_identical(stages$bootstrap$method, "inid")
    expect_identical(stages$bootstrap$B, 13L)
    expect_equal(stages$plot_data, candidate, info = case)
  }
})

test_that("npudens and npudist staged plot data can be rendered without estimator re-entry", {
  dens.proto <- getFromNamespace(".np_plot_proto_npudens_fixed_none_data", "np")
  dist.proto <- getFromNamespace(".np_plot_proto_npudist_fixed_none_data", "np")
  render <- getFromNamespace(".np_plot_proto_rectangular_surface_base_render", "np")
  withr::local_options(np.messages = FALSE)
  set.seed(2403)

  n <- 55L
  x <- data.frame(x1 = rnorm(n), x2 = runif(n, -1, 1))
  dens.bw <- npudensbw(
    dat = x,
    bws = c(0.55, 0.60),
    bandwidth.compute = FALSE,
    bwtype = "fixed"
  )
  dist.bw <- npudistbw(
    dat = x,
    bws = c(0.55, 0.60),
    bandwidth.compute = FALSE,
    bwtype = "fixed"
  )
  dens.data <- dens.proto(dens.bw, xdat = x, neval = 5L)
  dist.data <- dist.proto(dist.bw, xdat = x, neval = 5L)
  pdf.file <- tempfile(fileext = ".pdf")
  grDevices::pdf(pdf.file)
  on.exit(if (grDevices::dev.cur() > 1L) grDevices::dev.off(), add = TRUE)

  expect_identical(render(dens.data, perspective = FALSE), dens.data)
  expect_identical(render(dens.data, perspective = TRUE), dens.data)
  expect_identical(render(dist.data, perspective = FALSE), dist.data)
  expect_identical(render(dist.data, perspective = TRUE), dist.data)
  grDevices::dev.off()
  expect_true(file.exists(pdf.file))
  expect_gt(file.info(pdf.file)$size, 0)
})
