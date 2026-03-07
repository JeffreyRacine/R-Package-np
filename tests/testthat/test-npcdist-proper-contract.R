library(np)

test_that("proper helper projection enforces monotone bounded slices", {
  p.fun <- getFromNamespace(".np_condist_project_bounded_isotonic", "np")
  w.fun <- getFromNamespace(".np_condens_trapezoid_weights", "np")

  w <- w.fun(c(0, 1, 2, 3))
  f <- c(-0.2, 0.7, 0.4, 1.2)
  g <- p.fun(f = f, w = w, tol = 1e-12)

  expect_true(all(g >= -1e-10))
  expect_true(all(g <= 1 + 1e-10))
  expect_true(all(diff(g) >= -1e-10))
})

test_that("proper helpers support shadow-validation metadata path", {
  old.opt <- getOption("np.condist.proper.shadow")
  options(np.condist.proper.shadow = TRUE)
  on.exit(options(np.condist.proper.shadow = old.opt), add = TRUE)

  set.seed(20260307)
  x <- runif(30)
  y <- rnorm(30)
  y.grid <- seq(min(y) - 0.2, max(y) + 0.2, length.out = 15)
  x.grid <- c(0.2, 0.8)
  nd <- do.call(rbind, lapply(x.grid, function(xx) data.frame(y = y.grid, x = xx)))

  bw <- npcdistbw(
    xdat = data.frame(x = x),
    ydat = data.frame(y = y),
    bws = c(0.3, 0.3),
    bandwidth.compute = FALSE,
    regtype = "lp",
    degree = 3L
  )

  fit <- npcdist(
    bws = bw,
    txdat = data.frame(x = x),
    tydat = data.frame(y = y),
    exdat = nd["x"],
    eydat = nd["y"],
    proper = TRUE
  )

  expect_true(isTRUE(fit$proper.requested))
  expect_false(isTRUE(fit$proper.applied))
  expect_identical(fit$proper.info$reason, "repair_disabled_shadow_mode")
})

test_that("proper=FALSE preserves legacy npcdist output on explicit grids", {
  set.seed(20260307)
  x <- runif(50, -1, 1)
  y <- sin(2 * pi * x) + rnorm(50, sd = 0.2)
  y.grid <- seq(min(y) - 0.3, max(y) + 0.3, length.out = 25)
  x.grid <- c(-0.5, 0.5)
  nd <- do.call(rbind, lapply(x.grid, function(xx) data.frame(y = y.grid, x = xx)))

  bw <- npcdistbw(
    xdat = data.frame(x = x),
    ydat = data.frame(y = y),
    bws = c(0.28, 0.22),
    bandwidth.compute = FALSE,
    regtype = "lp",
    degree = 3L
  )

  fit.base <- npcdist(
    bws = bw,
    txdat = data.frame(x = x),
    tydat = data.frame(y = y),
    exdat = nd["x"],
    eydat = nd["y"]
  )
  fit.flag <- npcdist(
    bws = bw,
    txdat = data.frame(x = x),
    tydat = data.frame(y = y),
    exdat = nd["x"],
    eydat = nd["y"],
    proper = FALSE
  )

  expect_equal(fit.flag$condist, fit.base$condist, tolerance = 1e-12)
  expect_false(isTRUE(fit.flag$proper.requested))
  expect_false(isTRUE(fit.flag$proper.applied))
})

test_that("proper explicit-grid fit repairs slices and preserves raw values", {
  set.seed(2)
  n <- 80L
  x <- runif(n, -1, 1)
  y <- sin(2 * pi * x) + rnorm(n, sd = 0.20)
  x.grid <- unname(quantile(x, probs = c(0.2, 0.5, 0.8)))
  y.grid <- seq(min(y) - 0.3, max(y) + 0.3, length.out = 100L)
  nd <- do.call(rbind, lapply(x.grid, function(xx) data.frame(y = y.grid, x = xx)))

  bw <- npcdistbw(
    xdat = data.frame(x = x),
    ydat = data.frame(y = y),
    bws = c(0.28, 0.22),
    bandwidth.compute = FALSE,
    regtype = "lp",
    degree = 3L
  )

  fit.raw <- npcdist(
    bws = bw,
    txdat = data.frame(x = x),
    tydat = data.frame(y = y),
    exdat = nd["x"],
    eydat = nd["y"]
  )
  fit.proper <- npcdist(
    bws = bw,
    txdat = data.frame(x = x),
    tydat = data.frame(y = y),
    exdat = nd["x"],
    eydat = nd["y"],
    proper = TRUE
  )

  expect_true(isTRUE(fit.proper$proper.requested))
  expect_true(isTRUE(fit.proper$proper.applied))
  expect_identical(fit.proper$proper.method, "isotonic")
  expect_equal(fit.proper$condist.raw, fit.raw$condist, tolerance = 1e-12)
  expect_true(all(fit.proper$condist >= -1e-8))
  expect_true(all(fit.proper$condist <= 1 + 1e-8))

  split.idx <- split(seq_len(nrow(nd)), nd$x)
  for (idx in split.idx)
    expect_true(all(diff(fit.proper$condist[idx]) >= -1e-8))

  expect_error(se(fit.proper), "unavailable for repaired conditional distributions")
  expect_true(any(grepl("Proper distribution repair", capture.output(print(fit.proper)), fixed = TRUE)))
  expect_true(any(grepl("Proper distribution repair", capture.output(summary(fit.proper)), fixed = TRUE)))
})

test_that("proper request on paired evaluation stores metadata without altering fit", {
  set.seed(20260307)
  x <- runif(60)
  y <- rnorm(60)
  bw <- npcdistbw(
    xdat = data.frame(x = x),
    ydat = data.frame(y = y),
    bws = c(0.25, 0.25),
    bandwidth.compute = FALSE
  )

  fit.raw <- npcdist(
    bws = bw,
    txdat = data.frame(x = x),
    tydat = data.frame(y = y)
  )
  fit.req <- npcdist(
    bws = bw,
    txdat = data.frame(x = x),
    tydat = data.frame(y = y),
    proper = TRUE
  )

  expect_true(isTRUE(fit.req$proper.requested))
  expect_false(isTRUE(fit.req$proper.applied))
  expect_identical(fit.req$proper.info$reason, "no_eval_grid")
  expect_equal(fit.req$condist, fit.raw$condist, tolerance = 1e-12)
})

test_that("predict inherits proper request and validates newdata geometry", {
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
  fit.req <- npcdist(
    bws = bw,
    txdat = data.frame(x = x),
    tydat = data.frame(y = y),
    proper = TRUE
  )

  x.grid <- unname(quantile(x, probs = c(0.25, 0.75)))
  y.grid <- seq(min(y) - 0.2, max(y) + 0.2, length.out = 50L)
  nd.good <- do.call(rbind, lapply(x.grid, function(xx) data.frame(y = y.grid, x = xx)))
  pred <- predict(fit.req, newdata = nd.good)
  expect_true(all(pred >= -1e-8))
  expect_true(all(pred <= 1 + 1e-8))

  split.idx <- split(seq_len(nrow(nd.good)), nd.good$x)
  for (idx in split.idx)
    expect_true(all(diff(pred[idx]) >= -1e-8))

  nd.bad <- data.frame(y = y.grid[1:4], x = c(0.1, 0.2, 0.3, 0.4))
  expect_error(
    predict(fit.req, newdata = nd.bad),
    "requires repeated fixed-x slices|requires every fixed-x slice"
  )
})
