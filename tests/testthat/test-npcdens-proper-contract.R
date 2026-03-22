library(np)

test_that("proper helper weights and projection satisfy core invariants", {
  w.fun <- getFromNamespace(".np_condens_trapezoid_weights", "np")
  p.fun <- getFromNamespace(".np_condens_project_weighted_simplex", "np")

  w <- w.fun(c(0, 1, 2))
  expect_equal(w, c(0.5, 1, 0.5), tolerance = 1e-12)

  f <- c(-0.2, 0.4, 1.8)
  g <- p.fun(f = f, w = w, mass = 1, tol = 1e-12)
  expect_true(all(g >= -1e-10))
  expect_equal(sum(w * g), 1, tolerance = 1e-8)
})

test_that("proper helpers support shadow-validation metadata path", {
  old.opt <- getOption("np.condens.proper.shadow")
  options(np.condens.proper.shadow = TRUE)
  on.exit(options(np.condens.proper.shadow = old.opt), add = TRUE)

  set.seed(20260307)
  x <- runif(30)
  y <- rnorm(30)
  y.grid <- seq(min(y) - 0.2, max(y) + 0.2, length.out = 15)
  x.grid <- c(0.2, 0.8)
  nd <- do.call(rbind, lapply(x.grid, function(xx) data.frame(y = y.grid, x = xx)))

  bw <- npcdensbw(
    xdat = data.frame(x = x),
    ydat = data.frame(y = y),
    bws = c(0.3, 0.3),
    bandwidth.compute = FALSE,
    regtype = "lp",
    degree = 3L
  )

  fit <- npcdens(
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

test_that("proper=FALSE preserves legacy npcdens output on explicit grids", {
  set.seed(20260307)
  x <- runif(50, -1, 1)
  y <- sin(2 * pi * x) + rnorm(50, sd = 0.2)
  y.grid <- seq(min(y) - 0.3, max(y) + 0.3, length.out = 25)
  x.grid <- c(-0.5, 0.5)
  nd <- do.call(rbind, lapply(x.grid, function(xx) data.frame(y = y.grid, x = xx)))

  bw <- npcdensbw(
    xdat = data.frame(x = x),
    ydat = data.frame(y = y),
    bws = c(0.28, 0.22),
    bandwidth.compute = FALSE,
    regtype = "lp",
    degree = 3L
  )

  fit.base <- npcdens(
    bws = bw,
    txdat = data.frame(x = x),
    tydat = data.frame(y = y),
    exdat = nd["x"],
    eydat = nd["y"]
  )
  fit.flag <- npcdens(
    bws = bw,
    txdat = data.frame(x = x),
    tydat = data.frame(y = y),
    exdat = nd["x"],
    eydat = nd["y"],
    proper = FALSE
  )

  expect_equal(fit.flag$condens, fit.base$condens, tolerance = 1e-12)
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

  bw <- npcdensbw(
    xdat = data.frame(x = x),
    ydat = data.frame(y = y),
    bws = c(0.28, 0.22),
    bandwidth.compute = FALSE,
    regtype = "lp",
    degree = 3L
  )

  fit.raw <- npcdens(
    bws = bw,
    txdat = data.frame(x = x),
    tydat = data.frame(y = y),
    exdat = nd["x"],
    eydat = nd["y"]
  )
  fit.proper <- npcdens(
    bws = bw,
    txdat = data.frame(x = x),
    tydat = data.frame(y = y),
    exdat = nd["x"],
    eydat = nd["y"],
    proper = TRUE
  )

  expect_true(isTRUE(fit.proper$proper.requested))
  expect_true(isTRUE(fit.proper$proper.applied))
  expect_identical(fit.proper$proper.method, "project")
  expect_equal(fit.proper$condens.raw, fit.raw$condens, tolerance = 1e-12)
  expect_true(all(fit.proper$condens >= -1e-8))

  w.fun <- getFromNamespace(".np_condens_trapezoid_weights", "np")
  w <- w.fun(y.grid)
  split.idx <- split(seq_len(nrow(nd)), nd$x)
  mass <- vapply(split.idx, function(idx) sum(w * fit.proper$condens[idx]), numeric(1))
  expect_equal(unname(mass), rep(1, length(mass)), tolerance = 1e-6)

  expect_error(se(fit.proper), "unavailable for repaired conditional densities")
  expect_true(any(grepl("Proper density repair", capture.output(print(fit.proper)), fixed = TRUE)))
  expect_true(any(grepl("Proper density repair", capture.output(summary(fit.proper)), fixed = TRUE)))
})

test_that("proper request is a no-op for local-constant explicit-grid fits", {
  set.seed(20260322)
  x <- runif(60, -1, 1)
  y <- sin(2 * pi * x) + rnorm(60, sd = 0.15)
  y.grid <- seq(min(y) - 0.2, max(y) + 0.2, length.out = 40L)
  x.grid <- c(-0.4, 0.4)
  nd <- do.call(rbind, lapply(x.grid, function(xx) data.frame(y = y.grid, x = xx)))

  bw <- npcdensbw(
    xdat = data.frame(x = x),
    ydat = data.frame(y = y),
    bws = c(0.3, 0.24),
    bandwidth.compute = FALSE,
    regtype = "lc"
  )

  fit.raw <- npcdens(
    bws = bw,
    txdat = data.frame(x = x),
    tydat = data.frame(y = y),
    exdat = nd["x"],
    eydat = nd["y"]
  )
  fit.req <- npcdens(
    bws = bw,
    txdat = data.frame(x = x),
    tydat = data.frame(y = y),
    exdat = nd["x"],
    eydat = nd["y"],
    proper = TRUE
  )

  expect_true(isTRUE(fit.req$proper.requested))
  expect_false(isTRUE(fit.req$proper.applied))
  expect_identical(fit.req$proper.info$reason, "already_proper")
  expect_equal(fit.req$condens, fit.raw$condens, tolerance = 1e-12)
  expect_null(fit.req$condens.raw)
})

test_that("proper request is a no-op for degree-zero lp explicit-grid fits", {
  set.seed(20260322)
  x <- runif(60, -1, 1)
  y <- cos(2 * pi * x) + rnorm(60, sd = 0.15)
  y.grid <- seq(min(y) - 0.2, max(y) + 0.2, length.out = 40L)
  x.grid <- c(-0.25, 0.25)
  nd <- do.call(rbind, lapply(x.grid, function(xx) data.frame(y = y.grid, x = xx)))

  bw <- npcdensbw(
    xdat = data.frame(x = x),
    ydat = data.frame(y = y),
    bws = c(0.27, 0.21),
    bandwidth.compute = FALSE,
    regtype = "lp",
    degree = 0L
  )

  fit.raw <- npcdens(
    bws = bw,
    txdat = data.frame(x = x),
    tydat = data.frame(y = y),
    exdat = nd["x"],
    eydat = nd["y"]
  )
  fit.req <- npcdens(
    bws = bw,
    txdat = data.frame(x = x),
    tydat = data.frame(y = y),
    exdat = nd["x"],
    eydat = nd["y"],
    proper = TRUE
  )

  expect_true(isTRUE(fit.req$proper.requested))
  expect_false(isTRUE(fit.req$proper.applied))
  expect_identical(fit.req$proper.info$reason, "already_proper")
  expect_equal(fit.req$condens, fit.raw$condens, tolerance = 1e-12)
  expect_null(fit.req$condens.raw)
})

test_that("proper request on paired evaluation stores metadata without altering fit", {
  set.seed(20260307)
  x <- runif(60)
  y <- rnorm(60)
  bw <- npcdensbw(
    xdat = data.frame(x = x),
    ydat = data.frame(y = y),
    bws = c(0.25, 0.25),
    bandwidth.compute = FALSE,
    regtype = "lp",
    degree = 3L
  )

  fit.raw <- npcdens(
    bws = bw,
    txdat = data.frame(x = x),
    tydat = data.frame(y = y)
  )
  fit.req <- npcdens(
    bws = bw,
    txdat = data.frame(x = x),
    tydat = data.frame(y = y),
    proper = TRUE
  )

  expect_true(isTRUE(fit.req$proper.requested))
  expect_false(isTRUE(fit.req$proper.applied))
  expect_identical(fit.req$proper.info$reason, "no_eval_grid")
  expect_equal(fit.req$condens, fit.raw$condens, tolerance = 1e-12)
})

test_that("predict inherits proper request and validates newdata geometry", {
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
  fit.req <- npcdens(
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

  nd.bad <- data.frame(y = y.grid[1:4], x = c(0.1, 0.2, 0.3, 0.4))
  expect_error(
    predict(fit.req, newdata = nd.bad),
    "requires repeated fixed-x slices|requires every fixed-x slice"
  )
})
