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
    plot.behavior = "data",
    perspective = TRUE,
    view = "fixed"
  ))

  expect_type(out, "list")
  expect_true("cd1" %in% names(out))
  expect_s3_class(out$cd1, "condistribution")
  expect_true(isTRUE(out$cd1$proper.requested))
  expect_true(isTRUE(out$cd1$proper.applied))
  expect_true(all(out$cd1$condist >= -1e-8))
  expect_true(all(out$cd1$condist <= 1 + 1e-8))

  split.idx <- split(seq_along(out$cd1$condist), out$cd1$xeval)
  for (idx in split.idx)
    expect_true(all(diff(out$cd1$condist[idx]) >= -1e-8))
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
    plot.behavior = "data",
    perspective = FALSE
  ))

  expect_length(out, 2L)
  expect_s3_class(out[[1]], "condistribution")
  expect_s3_class(out[[2]], "condistribution")
  expect_true(isTRUE(out[[1]]$proper.requested))
  expect_false(isTRUE(out[[1]]$proper.applied))
  expect_true(isTRUE(out[[2]]$proper.requested))
  expect_true(isTRUE(out[[2]]$proper.applied))
  expect_true(all(out[[2]]$condist >= -1e-8))
  expect_true(all(out[[2]]$condist <= 1 + 1e-8))
  expect_true(all(diff(out[[2]]$condist) >= -1e-8))
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
    plot.behavior = "data",
    perspective = TRUE,
    view = "fixed",
    plot.errors.method = "asymptotic"
  ))
  boot <- suppressWarnings(plot(
    fit,
    plot.behavior = "data",
    perspective = TRUE,
    view = "fixed",
    plot.errors.method = "bootstrap",
    plot.errors.boot.method = "inid",
    plot.errors.boot.num = 19,
    plot.errors.type = "pointwise"
  ))

  expect_true(isTRUE(asym$cd1$proper.applied))
  expect_true(all(asym$cd1$condist >= -1e-8))
  expect_true(all(asym$cd1$condist <= 1 + 1e-8))
  expect_true(all(is.finite(asym$cd1$conderr)))

  expect_true(isTRUE(boot$cd1$proper.applied))
  expect_true(all(boot$cd1$condist >= -1e-8))
  expect_true(all(boot$cd1$condist <= 1 + 1e-8))
  expect_true(all(is.finite(boot$cd1$conderr)))

  xeval <- if (is.data.frame(boot$cd1$xeval)) boot$cd1$xeval[[1L]] else as.vector(boot$cd1$xeval)
  split.idx <- split(seq_along(boot$cd1$condist), xeval)
  for (idx in split.idx)
    expect_true(all(diff(boot$cd1$condist[idx]) >= -1e-8))
})

test_that("plot condistribution bootstrap supports mixed-bound fixed helper path", {
  set.seed(42)
  n <- 60L
  x <- runif(n)
  y <- runif(n)

  fit <- npcdist(y ~ x, cykerbound = "range")
  out <- suppressWarnings(plot(
    fit,
    plot.behavior = "data",
    view = "fixed",
    plot.data.overlay = FALSE,
    plot.errors.method = "bootstrap",
    plot.errors.boot.method = "inid",
    plot.errors.boot.num = 5L,
    plot.errors.type = "pointwise"
  ))

  expect_s3_class(out$cd1, "condistribution")
  expect_true(all(is.finite(out$cd1$condist)))
  expect_true(all(is.finite(out$cd1$conderr)))
})
