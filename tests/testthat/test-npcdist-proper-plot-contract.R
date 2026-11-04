library(np)

test_that("plot condistribution 2D data payload applies proper repair on supported grids", {
  set.seed(2)
  n <- 80L
  x <- runif(n, -1, 1)
  y <- sin(2 * pi * x) + rnorm(n, sd = 0.20)

  bw <- npcdistbw(
    xdat = data.frame(x = x),
    ydat = data.frame(y = y),
    bws = c(0.28, 0.22),
    bandwidth.compute = FALSE,
    regtype = "lp",
    degree = 3L
  )
  fit <- npcdist(
    bws = bw,
    txdat = data.frame(x = x),
    tydat = data.frame(y = y),
    proper = TRUE
  )

  out <- suppressWarnings(plot(
    fit,
    output = "data",
    perspective = TRUE,
    view = "fixed"
  ))

  expect_type(out, "list")
  expect_true("panel_1" %in% names(out$panels))
  expect_s3_class(out, "np_plot_data")
  expect_true(isTRUE(out$inference[[1L]]$projection$requested))
  expect_true(isTRUE(out$inference[[1L]]$projection$applied))
  expect_true(all(out$panels[[1L]]$values$estimate >= -1e-8))
  expect_true(all(out$panels[[1L]]$values$estimate <= 1 + 1e-8))

  split.idx <- split(seq_along(out$panels[[1L]]$values$estimate), plot_coordinates(out, 1L, "X"))
  for (idx in split.idx)
    expect_true(all(diff(out$panels[[1L]]$values$estimate[idx]) >= -1e-8))
})

test_that("plot condistribution 1D data payload repairs only y-varying panels", {
  set.seed(2)
  n <- 80L
  x <- runif(n, -1, 1)
  y <- sin(2 * pi * x) + rnorm(n, sd = 0.20)

  bw <- npcdistbw(
    xdat = data.frame(x = x),
    ydat = data.frame(y = y),
    bws = c(0.28, 0.22),
    bandwidth.compute = FALSE,
    regtype = "lp",
    degree = 3L
  )
  fit <- npcdist(
    bws = bw,
    txdat = data.frame(x = x),
    tydat = data.frame(y = y),
    proper = TRUE
  )

  out <- suppressWarnings(plot(
    fit,
    output = "data",
    perspective = FALSE
  ))

  expect_length(out$panels, 2L)
  expect_s3_class(out, "np_plot_data")
  expect_s3_class(out, "np_plot_data")
  expect_true(isTRUE(out$inference[[1L]]$projection$requested))
  expect_false(isTRUE(out$inference[[1L]]$projection$applied))
  expect_true(isTRUE(out$inference[[2L]]$projection$requested))
  expect_true(isTRUE(out$inference[[2L]]$projection$applied))
  expect_true(all(out$panels[[2L]]$values$estimate >= -1e-8))
  expect_true(all(out$panels[[2L]]$values$estimate <= 1 + 1e-8))
  expect_true(all(diff(out$panels[[2L]]$values$estimate) >= -1e-8))
})

test_that("plot condistribution supports asymptotic and bootstrap errors on proper grids", {
  set.seed(2)
  x <- runif(60, -1, 1)
  y <- sin(2 * pi * x) + rnorm(60, sd = 0.2)

  bw <- npcdistbw(
    xdat = data.frame(x = x),
    ydat = data.frame(y = y),
    bws = c(0.28, 0.22),
    bandwidth.compute = FALSE,
    regtype = "lp",
    degree = 3L
  )
  fit <- npcdist(
    bws = bw,
    txdat = data.frame(x = x),
    tydat = data.frame(y = y),
    proper = TRUE
  )

  asym <- suppressWarnings(plot(
    fit,
    output = "data",
    perspective = TRUE,
    view = "fixed",
    errors = "asymptotic"
  ))
  boot <- suppressWarnings(plot(
    fit,
    output = "data",
    perspective = TRUE,
    view = "fixed",
    errors = "bootstrap",
    bootstrap = "inid",
    B = 19,
    band = "pointwise"
  ))

  expect_true(isTRUE(asym$inference[[1L]]$projection$applied))
  expect_true(all(asym$panels[[1L]]$values$estimate >= -1e-8))
  expect_true(all(asym$panels[[1L]]$values$estimate <= 1 + 1e-8))
  expect_true(all(is.finite(as.matrix(asym$panels[[1L]]$intervals[c("lower", "upper")]))))

  expect_true(isTRUE(boot$inference[[1L]]$projection$applied))
  expect_true(all(boot$panels[[1L]]$values$estimate >= -1e-8))
  expect_true(all(boot$panels[[1L]]$values$estimate <= 1 + 1e-8))
  expect_true(all(is.finite(as.matrix(boot$panels[[1L]]$intervals[c("lower", "upper")]))))

  xeval <- if (is.data.frame(plot_coordinates(boot, 1L, "X"))) plot_coordinates(boot, 1L, "X")[[1L]] else as.vector(plot_coordinates(boot, 1L, "X"))
  split.idx <- split(seq_along(boot$panels[[1L]]$values$estimate), xeval)
  for (idx in split.idx)
    expect_true(all(diff(boot$panels[[1L]]$values$estimate[idx]) >= -1e-8))
})

test_that("plot condistribution bootstrap supports mixed-bound fixed helper path", {
  set.seed(42)
  n <- 60L
  x <- runif(n)
  y <- runif(n)

  fit <- npcdist(y ~ x, cykerbound = "range")
  out <- suppressWarnings(plot(
    fit,
    output = "data",
    view = "fixed",
    data_overlay = FALSE,
    errors = "bootstrap",
    bootstrap = "inid",
    B = 5L,
    band = "pointwise"
  ))

  expect_s3_class(out, "np_plot_data")
  expect_true(all(is.finite(out$panels[[1L]]$values$estimate)))
  expect_true(all(is.finite(as.matrix(out$panels[[1L]]$intervals[c("lower", "upper")]))))
})
