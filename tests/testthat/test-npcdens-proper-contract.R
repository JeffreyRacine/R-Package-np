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

test_that("slice mode can properize fitted npcdens values when apply='fitted'", {
  set.seed(20260322)
  x <- runif(36, -1, 1)
  y <- sin(2 * pi * x) + rnorm(36, sd = 0.18)

  bw <- npcdensbw(
    xdat = data.frame(x = x),
    ydat = data.frame(y = y),
    bws = c(0.29, 0.23),
    bandwidth.compute = FALSE,
    regtype = "lp",
    degree = 3L
  )

  ctrl <- list(mode = "slice", apply = "fitted", slice.grid.size = 31L, slice.extend.factor = 0)

  fit.raw <- npcdens(
    bws = bw,
    txdat = data.frame(x = x),
    tydat = data.frame(y = y)
  )
  fit.slice <- npcdens(
    bws = bw,
    txdat = data.frame(x = x),
    tydat = data.frame(y = y),
    proper = TRUE,
    proper.control = ctrl
  )

  build.grid <- getFromNamespace(".np_condens_build_slice_eval_grid", "np")
  grid.eval <- build.grid(
    object = fit.raw,
    slice.context = list(
      txdat = data.frame(x = x),
      tydat = data.frame(y = y),
      exdat = data.frame(x = x),
      eydat = data.frame(y = y)
    ),
    proper.control = ctrl
  )

  grid.fit <- npcdens(
    bws = bw,
    txdat = data.frame(x = x),
    tydat = data.frame(y = y),
    exdat = grid.eval$exdat,
    eydat = grid.eval$eydat
  )
  grid.proper <- npcdens(
    bws = bw,
    txdat = data.frame(x = x),
    tydat = data.frame(y = y),
    exdat = grid.eval$exdat,
    eydat = grid.eval$eydat,
    proper = TRUE
  )

  ypos <- match(y, grid.eval$y.grid)
  oracle <- numeric(length(y))
  for (i in seq_along(grid.eval$groups)) {
    idx.req <- grid.eval$groups[[i]]
    idx.grid <- grid.eval$grid.slices[[i]]
    oracle[idx.req] <- grid.proper$condens[idx.grid[ypos[idx.req]]]
  }

  expect_true(isTRUE(fit.slice$proper.requested))
  expect_true(isTRUE(fit.slice$proper.applied))
  expect_identical(fit.slice$proper.info$route, "slice")
  expect_identical(fit.slice$proper.info$apply.scope, "fitted")
  expect_equal(fit.slice$condens.raw, fit.raw$condens, tolerance = 1e-12)
  expect_equal(fit.slice$condens, oracle, tolerance = 1e-10)
  expect_true(any(abs(fit.slice$condens - fit.raw$condens) > 1e-8))
})

test_that("apply='fitted' leaves explicit-evaluation npcdens objects unchanged", {
  set.seed(20260322)
  x <- runif(70, -1, 1)
  y <- cos(2 * pi * x) + rnorm(70, sd = 0.2)
  nd <- rbind(
    data.frame(y = c(-0.8, -0.1, 0.55), x = rep(-0.4, 3L)),
    data.frame(y = c(-0.35, 0.15, 0.8, 1.1), x = rep(0.45, 4L))
  )

  bw <- npcdensbw(
    xdat = data.frame(x = x),
    ydat = data.frame(y = y),
    bws = c(0.26, 0.2),
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
  fit.req <- npcdens(
    bws = bw,
    txdat = data.frame(x = x),
    tydat = data.frame(y = y),
    exdat = nd["x"],
    eydat = nd["y"],
    proper = TRUE,
    proper.control = list(mode = "slice", apply = "fitted")
  )

  expect_true(isTRUE(fit.req$proper.requested))
  expect_false(isTRUE(fit.req$proper.applied))
  expect_identical(fit.req$proper.info$reason, "scope_not_selected")
  expect_true(isTRUE(fit.req$proper.info$supported))
  expect_equal(fit.req$condens, fit.raw$condens, tolerance = 1e-12)
})

test_that("apply='fitted' also leaves exact-grid npcdens evaluation objects unchanged", {
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
    degree = 3L
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
    proper = TRUE,
    proper.control = list(mode = "slice", apply = "fitted")
  )

  expect_true(isTRUE(fit.req$proper.requested))
  expect_false(isTRUE(fit.req$proper.applied))
  expect_identical(fit.req$proper.info$reason, "scope_not_selected")
  expect_true(isTRUE(fit.req$proper.info$supported))
  expect_equal(fit.req$condens, fit.raw$condens, tolerance = 1e-12)
})

test_that("slice mode defers to exact-grid repair when common grids are supplied", {
  set.seed(20260322)
  x <- runif(70, -1, 1)
  y <- sin(2 * pi * x) + rnorm(70, sd = 0.18)
  y.grid <- seq(min(y) - 0.25, max(y) + 0.25, length.out = 60L)
  x.grid <- c(-0.35, 0.35)
  nd <- do.call(rbind, lapply(x.grid, function(xx) data.frame(y = y.grid, x = xx)))

  bw <- npcdensbw(
    xdat = data.frame(x = x),
    ydat = data.frame(y = y),
    bws = c(0.28, 0.22),
    bandwidth.compute = FALSE,
    regtype = "lp",
    degree = 3L
  )

  fit.grid <- npcdens(
    bws = bw,
    txdat = data.frame(x = x),
    tydat = data.frame(y = y),
    exdat = nd["x"],
    eydat = nd["y"],
    proper = TRUE
  )
  fit.slice <- npcdens(
    bws = bw,
    txdat = data.frame(x = x),
    tydat = data.frame(y = y),
    exdat = nd["x"],
    eydat = nd["y"],
    proper = TRUE,
    proper.control = list(mode = "slice")
  )

  expect_true(isTRUE(fit.slice$proper.applied))
  expect_equal(fit.slice$condens, fit.grid$condens, tolerance = 1e-12)
})

test_that("slice mode on paired rows matches the internal explicit-grid oracle", {
  set.seed(20260322)
  x <- runif(80, -1, 1)
  y <- cos(2 * pi * x) + rnorm(80, sd = 0.2)
  nd <- rbind(
    data.frame(y = c(-0.8, -0.1, 0.55), x = rep(-0.4, 3L)),
    data.frame(y = c(-0.35, 0.15, 0.8, 1.1), x = rep(0.45, 4L))
  )

  bw <- npcdensbw(
    xdat = data.frame(x = x),
    ydat = data.frame(y = y),
    bws = c(0.26, 0.2),
    bandwidth.compute = FALSE,
    regtype = "lp",
    degree = 3L
  )

  ctrl <- list(mode = "slice", slice.grid.size = 21L, slice.extend.factor = 0)

  fit.slice <- npcdens(
    bws = bw,
    txdat = data.frame(x = x),
    tydat = data.frame(y = y),
    exdat = nd["x"],
    eydat = nd["y"],
    proper = TRUE,
    proper.control = ctrl
  )

  build.grid <- getFromNamespace(".np_condens_build_slice_eval_grid", "np")
  grid.eval <- build.grid(
    object = fit.slice,
    slice.context = list(
      txdat = data.frame(x = x),
      tydat = data.frame(y = y),
      exdat = nd["x"],
      eydat = nd["y"]
    ),
    proper.control = ctrl
  )
  oracle.fit <- npcdens(
    bws = bw,
    txdat = data.frame(x = x),
    tydat = data.frame(y = y),
    exdat = grid.eval$exdat,
    eydat = grid.eval$eydat,
    proper = TRUE
  )

  ypos <- match(nd$y, grid.eval$y.grid)
  oracle.values <- numeric(nrow(nd))
  for (i in seq_along(grid.eval$groups)) {
    idx.req <- grid.eval$groups[[i]]
    idx.grid <- grid.eval$grid.slices[[i]]
    oracle.values[idx.req] <- oracle.fit$condens[idx.grid[ypos[idx.req]]]
  }

  expect_true(isTRUE(fit.slice$proper.applied))
  expect_identical(fit.slice$proper.info$route, "slice")
  expect_true(all(fit.slice$condens >= -1e-8))
  expect_equal(fit.slice$condens, oracle.values, tolerance = 1e-10)
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
