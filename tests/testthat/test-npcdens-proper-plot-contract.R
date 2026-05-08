library(np)

test_that("plot condensity 2D data payload applies proper repair on supported grids", {
  set.seed(2)
  n <- 80L
  x <- runif(n, -1, 1)
  y <- sin(2 * pi * x) + rnorm(n, sd = 0.20)

  bw <- npcdensbw(
    xdat = data.frame(x = x),
    ydat = data.frame(y = y),
    bws = c(0.28, 0.22),
    bandwidth.compute = FALSE,
    regtype = "lp",
    degree = 3L
  )
  fit <- npcdens(
    bws = bw,
    txdat = data.frame(x = x),
    tydat = data.frame(y = y),
    proper = TRUE
  )

  out <- plot(
    fit,
    output = "data",
    perspective = TRUE,
    view = "fixed"
  )

  expect_type(out, "list")
  expect_true("cd1" %in% names(out))
  expect_s3_class(out$cd1, "condensity")
  expect_true(isTRUE(out$cd1$proper.requested))
  expect_true(isTRUE(out$cd1$proper.applied))
  expect_true(all(out$cd1$condens >= -1e-8))
})

test_that("plot condensity 1D data payload repairs only y-varying panels", {
  set.seed(2)
  n <- 80L
  x <- runif(n, -1, 1)
  y <- sin(2 * pi * x) + rnorm(n, sd = 0.20)

  bw <- npcdensbw(
    xdat = data.frame(x = x),
    ydat = data.frame(y = y),
    bws = c(0.28, 0.22),
    bandwidth.compute = FALSE,
    regtype = "lp",
    degree = 3L
  )
  fit <- npcdens(
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

  expect_length(out, 2L)
  expect_s3_class(out[[1]], "condensity")
  expect_s3_class(out[[2]], "condensity")
  expect_true(isTRUE(out[[1]]$proper.requested))
  expect_false(isTRUE(out[[1]]$proper.applied))
  expect_true(isTRUE(out[[2]]$proper.requested))
  expect_true(isTRUE(out[[2]]$proper.applied))
  expect_true(all(out[[2]]$condens >= -1e-8))
})

test_that("plot condensity supports asymptotic and bootstrap errors on proper grids", {
  set.seed(2)
  x <- runif(60, -1, 1)
  y <- sin(2 * pi * x) + rnorm(60, sd = 0.2)

  bw <- npcdensbw(
    xdat = data.frame(x = x),
    ydat = data.frame(y = y),
    bws = c(0.28, 0.22),
    bandwidth.compute = FALSE,
    regtype = "lp",
    degree = 3L
  )
  fit <- npcdens(
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

  expect_true(isTRUE(asym$cd1$proper.applied))
  expect_true(all(asym$cd1$condens >= -1e-8))
  expect_true(all(is.finite(asym$cd1$conderr)))

  expect_true(isTRUE(boot$cd1$proper.applied))
  expect_true(all(boot$cd1$condens >= -1e-8))
  expect_true(all(is.finite(boot$cd1$conderr)))

  yeval <- if (is.data.frame(boot$cd1$yeval)) boot$cd1$yeval[[1L]] else as.vector(boot$cd1$yeval)
  xeval <- if (is.data.frame(boot$cd1$xeval)) boot$cd1$xeval[[1L]] else as.vector(boot$cd1$xeval)
  w <- getFromNamespace(".np_condens_trapezoid_weights", "np")(sort(unique(yeval)))
  mass <- vapply(split(seq_along(boot$cd1$condens), xeval), function(idx) {
    sum(w * boot$cd1$condens[idx])
  }, numeric(1))
  expect_equal(unname(mass), rep(1, length(mass)), tolerance = 1e-8)
})

test_that("plot condensity bootstrap supports mixed-bound fixed helper path", {
  set.seed(42)
  n <- 60L
  x <- runif(n)
  y <- runif(n)

  fit <- npcdens(y ~ x, cykerbound = "range")
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

  expect_s3_class(out$cd1, "condensity")
  expect_true(all(is.finite(out$cd1$condens)))
  expect_true(all(is.finite(out$cd1$conderr)))
})
