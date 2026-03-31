test_that("npscoefhat reproduces npscoef fitted values and supports matrix RHS", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  old.auto <- getOption("npRmpi.autodispatch", FALSE)
  on.exit(options(npRmpi.autodispatch = old.auto), add = TRUE)
  on.exit(close_mpi_slaves(), add = TRUE)
  options(npRmpi.autodispatch = TRUE)

  set.seed(2468)
  n <- 90
  x <- runif(n)
  z <- runif(n)
  y <- (x^2) * z + 0.3 * x + rnorm(n, sd = 0.05)

  bw <- npscoefbw(
    xdat = x,
    zdat = z,
    ydat = y,
    bws = 0.15,
    bandwidth.compute = FALSE
  )

  tx <- data.frame(x = x)
  tz <- data.frame(z = z)
  ex <- data.frame(x = seq(min(x), max(x), length.out = 40))
  ez <- data.frame(z = seq(min(z), max(z), length.out = 40))

  fit.train <- npscoef(
    bws = bw,
    txdat = tx,
    tydat = y,
    tzdat = tz,
    iterate = FALSE
  )
  fit.eval <- npscoef(
    bws = bw,
    txdat = tx,
    tydat = y,
    tzdat = tz,
    exdat = ex,
    ezdat = ez,
    iterate = FALSE
  )

  H.train <- npscoefhat(
    bws = bw,
    txdat = tx,
    tzdat = tz,
    output = "matrix",
    iterate = FALSE
  )
  H.eval <- npscoefhat(
    bws = bw,
    txdat = tx,
    tzdat = tz,
    exdat = ex,
    ezdat = ez,
    output = "matrix",
    iterate = FALSE
  )

  expect_equal(as.vector(H.train %*% y), as.vector(fit.train$mean), tolerance = 1e-8)
  expect_equal(as.vector(H.eval %*% y), as.vector(fit.eval$mean), tolerance = 1e-8)

  ystar <- cbind(y, y + 0.1)
  hy <- npscoefhat(
    bws = bw,
    txdat = tx,
    tzdat = tz,
    exdat = ex,
    ezdat = ez,
    y = ystar,
    output = "apply",
    iterate = FALSE
  )
  expect_true(isTRUE(all.equal(
    hy,
    H.eval %*% ystar,
    tolerance = 1e-10,
    check.attributes = FALSE
  )))
  H.loo.lc <- npscoefhat(
    bws = bw,
    txdat = tx,
    tzdat = tz,
    output = "matrix",
    iterate = FALSE,
    leave.one.out = TRUE
  )
  expect_gt(max(abs(H.train - H.loo.lc)), 1e-8)

  bw.ll <- npscoefbw(
    xdat = x, zdat = z, ydat = y,
    bws = 0.16, regtype = "ll", bandwidth.compute = FALSE
  )
  H0.ll <- npscoefhat(
    bws = bw.ll, txdat = tx, tzdat = tz,
    output = "matrix", iterate = FALSE, leave.one.out = FALSE
  )
  H1.ll <- npscoefhat(
    bws = bw.ll, txdat = tx, tzdat = tz,
    output = "matrix", iterate = FALSE, leave.one.out = TRUE
  )
  expect_gt(max(abs(H0.ll - H1.ll)), 1e-8)
  fit0.ll <- npscoef(
    bws = bw.ll, txdat = tx, tydat = y, tzdat = tz,
    iterate = FALSE, errors = FALSE, leave.one.out = FALSE
  )
  fit1.ll <- npscoef(
    bws = bw.ll, txdat = tx, tydat = y, tzdat = tz,
    iterate = FALSE, errors = FALSE, leave.one.out = TRUE
  )
  expect_gt(max(abs(fit0.ll$mean - fit1.ll$mean)), 1e-8)

  ex2 <- data.frame(x = seq(min(x), max(x), length.out = 20))
  ez2 <- data.frame(z = seq(min(z), max(z), length.out = 20))
  expect_error(
    npscoefhat(
      bws = bw.ll,
      txdat = tx, tzdat = tz,
      exdat = ex2, ezdat = ez2,
      output = "matrix",
      iterate = FALSE,
      leave.one.out = TRUE
    ),
    "requires evaluation 'z' data to match training 'z' data"
  )
})

test_that("npplreghat reproduces npplreg fitted values and supports matrix RHS", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  old.auto <- getOption("npRmpi.autodispatch", FALSE)
  on.exit(options(npRmpi.autodispatch = old.auto), add = TRUE)
  on.exit(close_mpi_slaves(), add = TRUE)
  options(npRmpi.autodispatch = TRUE)

  set.seed(97531)
  n <- 120
  x <- runif(n)
  z <- runif(n)
  y <- sin(z) + 2.0 * x + rnorm(n, sd = 0.05)

  bw <- npplregbw(
    xdat = x,
    zdat = z,
    ydat = y,
    bws = matrix(c(0.2, 0.2), nrow = 2, ncol = 1),
    bandwidth.compute = FALSE
  )

  tx <- data.frame(x = x)
  tz <- data.frame(z = z)
  ex <- data.frame(x = seq(min(x), max(x), length.out = 35))
  ez <- data.frame(z = seq(min(z), max(z), length.out = 35))

  fit.train <- npplreg(
    bws = bw,
    txdat = tx,
    tydat = y,
    tzdat = tz
  )
  fit.eval <- npplreg(
    bws = bw,
    txdat = tx,
    tydat = y,
    tzdat = tz,
    exdat = ex,
    ezdat = ez
  )

  H.train <- npplreghat(
    bws = bw,
    txdat = tx,
    tzdat = tz,
    output = "matrix"
  )
  H.eval <- npplreghat(
    bws = bw,
    txdat = tx,
    tzdat = tz,
    exdat = ex,
    ezdat = ez,
    output = "matrix"
  )

  expect_equal(as.vector(H.train %*% y), as.vector(fit.train$mean), tolerance = 1e-8)
  expect_equal(as.vector(H.eval %*% y), as.vector(fit.eval$mean), tolerance = 1e-8)

  ystar <- cbind(y, y + 0.05)
  hy <- npplreghat(
    bws = bw,
    txdat = tx,
    tzdat = tz,
    exdat = ex,
    ezdat = ez,
    y = ystar,
    output = "apply"
  )
  expect_equal(hy, H.eval %*% ystar, tolerance = 1e-10)
})

test_that("npplreghat generalized-nn apply matches npplreg means in session mode", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  old.auto <- getOption("npRmpi.autodispatch", FALSE)
  on.exit(options(npRmpi.autodispatch = old.auto), add = TRUE)
  on.exit(close_mpi_slaves(), add = TRUE)
  options(npRmpi.autodispatch = TRUE)

  set.seed(20260308)
  n <- 80
  x <- runif(n)
  z <- runif(n)
  y <- sin(z) + 2.0 * x + rnorm(n, sd = 0.05)

  tx <- data.frame(x = x)
  tz <- data.frame(z = z)
  ex <- data.frame(x = seq(min(x), max(x), length.out = 25))
  ez <- data.frame(z = seq(min(z), max(z), length.out = 25))

  bw <- npplregbw(
    xdat = tx,
    zdat = tz,
    ydat = y,
    regtype = "ll",
    bwtype = "generalized_nn",
    bws = matrix(c(2, 9), nrow = 2, ncol = 1),
    bandwidth.compute = FALSE
  )

  fit <- npplreg(
    bws = bw,
    txdat = tx,
    tydat = y,
    tzdat = tz,
    exdat = ex,
    ezdat = ez
  )
  hy <- npplreghat(
    bws = bw,
    txdat = tx,
    tzdat = tz,
    exdat = ex,
    ezdat = ez,
    y = y,
    output = "apply"
  )

  expect_equal(as.vector(hy), as.vector(fit$mean), tolerance = 1e-10)
})

test_that("npplreghat adaptive-nn apply matches npplreg means in session mode", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  old.auto <- getOption("npRmpi.autodispatch", FALSE)
  on.exit(options(npRmpi.autodispatch = old.auto), add = TRUE)
  on.exit(close_mpi_slaves(), add = TRUE)
  options(npRmpi.autodispatch = TRUE)

  set.seed(20260311)
  n <- 80
  x <- runif(n)
  z <- runif(n)
  y <- sin(z) + 2.0 * x + rnorm(n, sd = 0.05)

  tx <- data.frame(x = x)
  tz <- data.frame(z = z)
  ex <- data.frame(x = seq(min(x), max(x), length.out = 25))
  ez <- data.frame(z = seq(min(z), max(z), length.out = 25))

  bw <- npplregbw(
    xdat = tx,
    zdat = tz,
    ydat = y,
    regtype = "ll",
    bwtype = "adaptive_nn",
    bws = matrix(c(9, 9), nrow = 2, ncol = 1),
    bandwidth.compute = FALSE
  )

  fit <- npplreg(
    bws = bw,
    txdat = tx,
    tydat = y,
    tzdat = tz,
    exdat = ex,
    ezdat = ez
  )
  hy <- npplreghat(
    bws = bw,
    txdat = tx,
    tzdat = tz,
    exdat = ex,
    ezdat = ez,
    y = y,
    output = "apply"
  )

  expect_equal(as.vector(hy), as.vector(fit$mean), tolerance = 1e-10)
})

test_that("npplreg generalized-nn inid plot helper completes in session mode", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  old.auto <- getOption("npRmpi.autodispatch", FALSE)
  on.exit(options(npRmpi.autodispatch = old.auto), add = TRUE)
  on.exit(close_mpi_slaves(), add = TRUE)
  options(npRmpi.autodispatch = TRUE)

  set.seed(20260308)
  n <- 80
  x <- runif(n)
  z <- runif(n)
  y <- sin(z) + 2.0 * x + rnorm(n, sd = 0.05)

  tx <- data.frame(x = x)
  tz <- data.frame(z = z)

  bw <- npplregbw(
    xdat = tx,
    zdat = tz,
    ydat = y,
    regtype = "ll",
    bwtype = "generalized_nn",
    bws = matrix(c(2, 9), nrow = 2, ncol = 1),
    bandwidth.compute = FALSE
  )

  set.seed(20260308)
  out <- plot(
    bw,
    xdat = tx,
    ydat = y,
    zdat = tz,
    plot.behavior = "data",
    perspective = FALSE,
    plot.errors.method = "bootstrap",
    plot.errors.boot.method = "inid",
    plot.errors.boot.num = 5
  )

  expect_true(is.list(out))
  expect_equal(length(out), 2L)
  expect_true(all(is.finite(out[[1L]]$mean)))
  expect_true(all(is.finite(out[[1L]]$merr)))
  expect_true(all(is.finite(out[[2L]]$mean)))
  expect_true(all(is.finite(out[[2L]]$merr)))
})

test_that("npplreg generalized-nn wild plot helper preserves means in session mode", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  old.auto <- getOption("npRmpi.autodispatch", FALSE)
  on.exit(options(npRmpi.autodispatch = old.auto), add = TRUE)
  on.exit(close_mpi_slaves(), add = TRUE)
  options(npRmpi.autodispatch = TRUE)

  set.seed(20260308)
  n <- 80
  x <- runif(n)
  z <- runif(n)
  y <- sin(z) + 2.0 * x + rnorm(n, sd = 0.05)

  tx <- data.frame(x = x)
  tz <- data.frame(z = z)

  bw <- npplregbw(
    xdat = tx,
    zdat = tz,
    ydat = y,
    regtype = "ll",
    bwtype = "generalized_nn",
    bws = matrix(c(2, 9), nrow = 2, ncol = 1),
    bandwidth.compute = FALSE
  )

  none.out <- plot(
    bw,
    xdat = tx,
    ydat = y,
    zdat = tz,
    plot.behavior = "data",
    perspective = FALSE,
    plot.errors.method = "none"
  )

  set.seed(20260308)
  wild.out <- plot(
    bw,
    xdat = tx,
    ydat = y,
    zdat = tz,
    plot.behavior = "data",
    perspective = FALSE,
    plot.errors.method = "bootstrap",
    plot.errors.boot.method = "wild",
    plot.errors.boot.num = 5
  )

  expect_equal(as.vector(wild.out[[1L]]$mean), as.vector(none.out[[1L]]$mean), tolerance = 1e-10)
  expect_equal(as.vector(wild.out[[2L]]$mean), as.vector(none.out[[2L]]$mean), tolerance = 1e-10)
  expect_true(all(is.finite(wild.out[[1L]]$merr)))
  expect_true(all(is.finite(wild.out[[2L]]$merr)))
})

test_that("npplreg adaptive-nn wild plot helper preserves means in session mode", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  old.auto <- getOption("npRmpi.autodispatch", FALSE)
  on.exit(options(npRmpi.autodispatch = old.auto), add = TRUE)
  on.exit(close_mpi_slaves(), add = TRUE)
  options(npRmpi.autodispatch = TRUE)

  set.seed(20260311)
  n <- 80
  x <- runif(n)
  z <- runif(n)
  y <- sin(z) + 2.0 * x + rnorm(n, sd = 0.05)

  tx <- data.frame(x = x)
  tz <- data.frame(z = z)

  bw <- npplregbw(
    xdat = tx,
    zdat = tz,
    ydat = y,
    regtype = "ll",
    bwtype = "adaptive_nn",
    bws = matrix(c(9, 9), nrow = 2, ncol = 1),
    bandwidth.compute = FALSE
  )

  none.out <- plot(
    bw,
    xdat = tx,
    ydat = y,
    zdat = tz,
    plot.behavior = "data",
    perspective = FALSE,
    plot.errors.method = "none"
  )

  set.seed(20260311)
  wild.out <- plot(
    bw,
    xdat = tx,
    ydat = y,
    zdat = tz,
    plot.behavior = "data",
    perspective = FALSE,
    plot.errors.method = "bootstrap",
    plot.errors.boot.method = "wild",
    plot.errors.boot.num = 5
  )

  expect_equal(as.vector(wild.out[[1L]]$mean), as.vector(none.out[[1L]]$mean), tolerance = 1e-10)
  expect_equal(as.vector(wild.out[[2L]]$mean), as.vector(none.out[[2L]]$mean), tolerance = 1e-10)
  expect_true(all(is.finite(wild.out[[1L]]$merr)))
  expect_true(all(is.finite(wild.out[[2L]]$merr)))
})

test_that("npplreg generalized-nn plot means match public estimator in session mode", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  old.auto <- getOption("npRmpi.autodispatch", FALSE)
  on.exit(options(npRmpi.autodispatch = old.auto), add = TRUE)
  on.exit(close_mpi_slaves(), add = TRUE)
  options(npRmpi.autodispatch = TRUE)

  set.seed(20260308)
  n <- 80
  x <- runif(n)
  z <- runif(n)
  y <- sin(z) + 2.0 * x + rnorm(n, sd = 0.05)

  tx <- data.frame(x = x)
  tz <- data.frame(z = z)

  bw <- npplregbw(
    xdat = tx,
    zdat = tz,
    ydat = y,
    regtype = "ll",
    bwtype = "generalized_nn",
    bws = matrix(c(2, 9), nrow = 2, ncol = 1),
    bandwidth.compute = FALSE
  )

  out <- plot(
    bw,
    xdat = tx,
    ydat = y,
    zdat = tz,
    plot.behavior = "data",
    perspective = FALSE,
    plot.errors.method = "none"
  )

  ex.x <- out[[1L]]$evalx
  ez.x <- out[[1L]]$evalz
  fit.x <- npplreg(
    bws = bw,
    txdat = tx,
    tydat = y,
    tzdat = tz,
    exdat = ex.x,
    ezdat = ez.x
  )

  ex.z <- out[[2L]]$evalx
  ez.z <- out[[2L]]$evalz
  fit.z <- npplreg(
    bws = bw,
    txdat = tx,
    tydat = y,
    tzdat = tz,
    exdat = ex.z,
    ezdat = ez.z
  )

  expect_equal(as.vector(out[[1L]]$mean), as.vector(fit.x$mean), tolerance = 1e-10)
  expect_equal(as.vector(out[[2L]]$mean), as.vector(fit.z$mean), tolerance = 1e-10)
})

test_that("npplreghat supports ll/lp with lp basis variants", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  old.auto <- getOption("npRmpi.autodispatch", FALSE)
  on.exit(options(npRmpi.autodispatch = old.auto), add = TRUE)
  on.exit(close_mpi_slaves(), add = TRUE)
  options(npRmpi.autodispatch = TRUE)

  set.seed(97532)
  n <- 90
  x <- runif(n)
  z <- runif(n)
  y <- sin(2 * z) + 1.25 * x + rnorm(n, sd = 0.04)
  tx <- data.frame(x = x)
  tz <- data.frame(z = z)
  ex <- data.frame(x = seq(min(x), max(x), length.out = 27))
  ez <- data.frame(z = seq(min(z), max(z), length.out = 27))

  cfgs <- list(
    list(regtype = "ll", basis = NULL, label = "ll"),
    list(regtype = "lp", basis = "glp", label = "lp-glp"),
    list(regtype = "lp", basis = "additive", label = "lp-additive"),
    list(regtype = "lp", basis = "tensor", label = "lp-tensor")
  )

  for (cfg in cfgs) {
    bw.args <- list(
      xdat = x,
      zdat = z,
      ydat = y,
      bws = matrix(c(0.22, 0.22), nrow = 2, ncol = 1),
      bandwidth.compute = FALSE,
      regtype = cfg$regtype
    )
    if (!is.null(cfg$basis)) {
      bw.args$basis <- cfg$basis
      bw.args$degree <- 2L
    }
    bw <- do.call(npplregbw, bw.args)
    if (!is.null(cfg$basis))
      expect_identical(bw$basis, cfg$basis, info = cfg$label)

    fit.eval <- npplreg(
      bws = bw,
      txdat = tx,
      tydat = y,
      tzdat = tz,
      exdat = ex,
      ezdat = ez
    )
    H.eval <- npplreghat(
      bws = bw,
      txdat = tx,
      tzdat = tz,
      exdat = ex,
      ezdat = ez,
      output = "matrix"
    )
    expect_equal(
      as.vector(H.eval %*% y),
      as.vector(fit.eval$mean),
      tolerance = 1e-8,
      info = cfg$label
    )
  }
})

test_that("npindexhat reproduces npindex fit and approximates gradient", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  old.auto <- getOption("npRmpi.autodispatch", FALSE)
  on.exit(options(npRmpi.autodispatch = old.auto), add = TRUE)
  on.exit(close_mpi_slaves(), add = TRUE)
  options(npRmpi.autodispatch = TRUE)

  set.seed(314159)
  n <- 110
  x1 <- runif(n)
  x2 <- runif(n)
  y <- sin(2 * (x1 + x2)) + rnorm(n, sd = 0.05)

  tx <- data.frame(x1 = x1, x2 = x2)
  bw <- npindexbw(xdat = tx, ydat = y, method = "ichimura", nmulti = 1)

  fit.mean <- npindex(
    bws = bw,
    txdat = tx,
    tydat = y,
    exdat = tx,
    gradients = FALSE
  )
  fit.grad <- npindex(
    bws = bw,
    txdat = tx,
    tydat = y,
    exdat = tx,
    gradients = TRUE
  )

  H0 <- npindexhat(
    bws = bw,
    txdat = tx,
    exdat = tx,
    s = 0L,
    output = "matrix"
  )
  H1 <- npindexhat(
    bws = bw,
    txdat = tx,
    exdat = tx,
    s = 1L,
    output = "matrix"
  )

  expect_equal(as.vector(H0 %*% y), as.vector(fit.mean$mean), tolerance = 1e-8)
  expect_equal(as.vector(H1 %*% y), as.vector(fit.grad$grad[, 1]), tolerance = 5e-3)

  ystar <- cbind(y, y - 0.05)
  hy <- npindexhat(
    bws = bw,
    txdat = tx,
    exdat = tx,
    y = ystar,
    s = 0L,
    output = "apply"
  )
  expect_true(isTRUE(all.equal(
    hy,
    H0 %*% ystar,
    tolerance = 1e-10,
    check.attributes = FALSE
  )))
})

test_that("npindex and npindexhat support ll/lp basis variants", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  old.auto <- getOption("npRmpi.autodispatch", FALSE)
  on.exit(options(npRmpi.autodispatch = old.auto), add = TRUE)
  on.exit(close_mpi_slaves(), add = TRUE)
  options(npRmpi.autodispatch = TRUE)

  set.seed(314160)
  n <- 90
  x1 <- runif(n)
  x2 <- runif(n)
  y <- sin(2 * (x1 + x2)) + rnorm(n, sd = 0.05)
  tx <- data.frame(x1 = x1, x2 = x2)
  ex <- tx[seq_len(30), , drop = FALSE]

  cfgs <- list(
    list(regtype = "ll", basis = NULL, label = "ll"),
    list(regtype = "lp", basis = "glp", label = "lp-glp"),
    list(regtype = "lp", basis = "additive", label = "lp-additive"),
    list(regtype = "lp", basis = "tensor", label = "lp-tensor")
  )

  for (cfg in cfgs) {
    bw.args <- list(
      xdat = tx,
      ydat = y,
      bws = c(1, 1, 0.25),
      bandwidth.compute = FALSE,
      regtype = cfg$regtype
    )
    if (!is.null(cfg$basis)) {
      bw.args$basis <- cfg$basis
      bw.args$degree <- 2L
    }
    bw <- do.call(npindexbw, bw.args)
    expect_identical(bw$regtype, cfg$regtype, info = cfg$label)
    if (!is.null(cfg$basis))
      expect_identical(bw$basis, cfg$basis, info = cfg$label)

    fit <- npindex(
      bws = bw,
      txdat = tx,
      tydat = y,
      exdat = ex,
      gradients = FALSE
    )
    H <- npindexhat(
      bws = bw,
      txdat = tx,
      exdat = ex,
      s = 0L,
      output = "matrix"
    )
    expect_equal(as.vector(H %*% y), as.vector(fit$mean), tolerance = 1e-8, info = cfg$label)
  }
})

test_that("npindexhat s=1 generalized-nn helper apply matches helper matrix in session mode", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  old.auto <- getOption("npRmpi.autodispatch", FALSE)
  on.exit(options(npRmpi.autodispatch = old.auto), add = TRUE)
  on.exit(close_mpi_slaves(), add = TRUE)
  options(npRmpi.autodispatch = TRUE)

  set.seed(20260308)
  n <- 70
  x1 <- runif(n)
  x2 <- runif(n)
  z <- runif(n)
  y <- (0.4 + x1) * sin(2 * pi * z) + 0.3 * x2 + rnorm(n, sd = 0.04)

  tx <- data.frame(x1 = x1, x2 = x2)
  ex <- tx[seq_len(18), , drop = FALSE]
  bw <- npindexbw(
    xdat = tx,
    ydat = y,
    regtype = "ll",
    bwtype = "generalized_nn",
    bandwidth.compute = FALSE,
    bws = c(1, 1, 9)
  )

  H1 <- npindexhat(
    bws = bw,
    txdat = tx,
    exdat = ex,
    s = 1L,
    output = "matrix"
  )
  hy <- npindexhat(
    bws = bw,
    txdat = tx,
    exdat = ex,
    y = y,
    s = 1L,
    output = "apply"
  )

  expect_equal(as.vector(hy), as.vector(H1 %*% y), tolerance = 1e-10)
})

test_that("npindex generalized-nn ll gradients match npindexhat s=1 helper in session mode", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  old.auto <- getOption("npRmpi.autodispatch", FALSE)
  on.exit(options(npRmpi.autodispatch = old.auto), add = TRUE)
  on.exit(close_mpi_slaves(), add = TRUE)
  options(npRmpi.autodispatch = TRUE)

  set.seed(20260308)
  n <- 70
  x1 <- runif(n)
  x2 <- runif(n)
  z <- runif(n)
  y <- (0.4 + x1) * sin(2 * pi * z) + 0.3 * x2 + rnorm(n, sd = 0.04)

  tx <- data.frame(x1 = x1, x2 = x2)
  ex <- tx[seq_len(18), , drop = FALSE]
  bw <- npindexbw(
    xdat = tx,
    ydat = y,
    regtype = "ll",
    bwtype = "generalized_nn",
    bandwidth.compute = FALSE,
    bws = c(1, 1, 9)
  )

  fit <- npindex(
    bws = bw,
    txdat = tx,
    tydat = y,
    exdat = ex,
    gradients = TRUE
  )
  hy <- npindexhat(
    bws = bw,
    txdat = tx,
    exdat = ex,
    y = y,
    s = 1L,
    output = "apply"
  )

  expect_equal(as.vector(fit$grad[, 1L]), as.vector(hy), tolerance = 1e-10)
})

test_that("npindexhat mean and derivative operators match core fits across bwtypes", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  old.auto <- getOption("npRmpi.autodispatch", FALSE)
  on.exit(options(npRmpi.autodispatch = old.auto), add = TRUE)
  on.exit(close_mpi_slaves(), add = TRUE)
  options(npRmpi.autodispatch = TRUE)

  set.seed(20260308)
  n <- 70
  x1 <- runif(n)
  x2 <- runif(n)
  y <- sin(x1 + x2) + rnorm(n, sd = 0.06)
  tx <- data.frame(x1 = x1, x2 = x2)
  ex <- tx[seq_len(20), , drop = FALSE]

  cfgs <- list(
    list(regtype = "lc", basis = NULL, degree = NULL, h = 0.85),
    list(regtype = "ll", basis = NULL, degree = NULL, h = 0.85),
    list(regtype = "lp", basis = "tensor", degree = 2L, h = 0.85)
  )

  for (bt in c("fixed", "generalized_nn", "adaptive_nn")) {
    for (cfg in cfgs) {
      bw.args <- list(
        xdat = tx,
        ydat = y,
        bws = c(1, 1, cfg$h),
        bandwidth.compute = FALSE,
        regtype = cfg$regtype,
        bwtype = bt
      )
      if (!is.null(cfg$basis)) {
        bw.args$basis <- cfg$basis
        bw.args$degree <- cfg$degree
      }
      bw <- do.call(npindexbw, bw.args)

      fit.mean <- npindex(
        bws = bw,
        txdat = tx,
        tydat = y,
        exdat = ex,
        gradients = FALSE
      )
      fit.grad <- npindex(
        bws = bw,
        txdat = tx,
        tydat = y,
        exdat = ex,
        gradients = TRUE
      )
      H.mean <- npindexhat(
        bws = bw,
        txdat = tx,
        exdat = ex,
        output = "matrix",
        s = 0L
      )
      a.mean <- npindexhat(
        bws = bw,
        txdat = tx,
        exdat = ex,
        y = y,
        output = "apply",
        s = 0L
      )
      H.grad <- npindexhat(
        bws = bw,
        txdat = tx,
        exdat = ex,
        output = "matrix",
        s = 1L
      )
      a.grad <- npindexhat(
        bws = bw,
        txdat = tx,
        exdat = ex,
        y = y,
        output = "apply",
        s = 1L
      )

      expect_equal(as.vector(H.mean %*% y), as.vector(fit.mean$mean), tolerance = 1e-8, info = paste(bt, cfg$regtype, "mean matrix"))
      expect_equal(as.vector(a.mean), as.vector(fit.mean$mean), tolerance = 1e-8, info = paste(bt, cfg$regtype, "mean apply"))
      expect_equal(as.vector(H.mean %*% y), as.vector(a.mean), tolerance = 1e-10, info = paste(bt, cfg$regtype, "mean matrix/apply"))
      expect_equal(as.vector(H.grad %*% y), as.vector(fit.grad$grad[, 1L]), tolerance = 1e-8, info = paste(bt, cfg$regtype, "grad matrix"))
      expect_equal(as.vector(a.grad), as.vector(fit.grad$grad[, 1L]), tolerance = 1e-8, info = paste(bt, cfg$regtype, "grad apply"))
      expect_equal(as.vector(H.grad %*% y), as.vector(a.grad), tolerance = 1e-10, info = paste(bt, cfg$regtype, "grad matrix/apply"))
    }
  }
})

test_that("npindex and npindexhat preserve nearest-neighbor bwtype semantics", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  old.auto <- getOption("npRmpi.autodispatch", FALSE)
  on.exit(options(npRmpi.autodispatch = old.auto), add = TRUE)
  on.exit(close_mpi_slaves(), add = TRUE)
  options(npRmpi.autodispatch = TRUE)

  set.seed(314161)
  n <- 70
  x1 <- runif(n)
  x2 <- runif(n)
  y <- sin(x1 + x2) + rnorm(n, sd = 0.06)
  tx <- data.frame(x1 = x1, x2 = x2)
  ex <- tx[seq_len(20), , drop = FALSE]

  cfgs <- list(
    list(regtype = "lc", basis = NULL, degree = NULL, h = 0.85),
    list(regtype = "ll", basis = NULL, degree = NULL, h = 0.85),
    list(regtype = "lp", basis = "tensor", degree = 2L, h = 0.85)
  )

  for (bt in c("generalized_nn", "adaptive_nn")) {
    for (cfg in cfgs) {
      bw.args <- list(
        xdat = tx,
        ydat = y,
        bws = c(1, 1, cfg$h),
        bandwidth.compute = FALSE,
        regtype = cfg$regtype,
        bwtype = bt
      )
      if (!is.null(cfg$basis)) {
        bw.args$basis <- cfg$basis
        bw.args$degree <- cfg$degree
      }
      bw <- do.call(npindexbw, bw.args)

      fit <- npindex(
        bws = bw,
        txdat = tx,
        tydat = y,
        exdat = ex,
        gradients = FALSE
      )
      H <- npindexhat(
        bws = bw,
        txdat = tx,
        exdat = ex,
        output = "matrix",
        s = 0L
      )
      hy <- npindexhat(
        bws = bw,
        txdat = tx,
        exdat = ex,
        y = y,
        output = "apply",
        s = 0L
      )

      expect_equal(as.vector(H %*% y), as.vector(fit$mean), tolerance = 1e-8, info = paste(bt, cfg$regtype, "matrix"))
      expect_equal(as.vector(hy), as.vector(fit$mean), tolerance = 1e-8, info = paste(bt, cfg$regtype, "apply"))
    }
  }

  sibw.bad.bounds <- npindexbw(
    xdat = tx,
    ydat = y,
    bws = c(1, 1, 0.22),
    bandwidth.compute = FALSE,
    regtype = "lp",
    basis = "tensor",
    degree = 2L,
    bernstein.basis = TRUE,
    bwtype = "fixed",
    ckertype = "epanechnikov",
    ckerorder = 2L,
    ckerbound = "fixed",
    ckerlb = 0.0,
    ckerub = 1.0
  )
  expect_error(
    npindexhat(
      bws = sibw.bad.bounds,
      txdat = tx,
      exdat = ex,
      s = 0L,
      output = "matrix"
    ),
    "Invalid bounds for 'ckerbound'|Evaluation data violate 'ckerbound' bounds"
  )
})

test_that("npscoefhat apply mode matches green core fits across bwtypes", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  old.auto <- getOption("npRmpi.autodispatch", FALSE)
  on.exit(options(npRmpi.autodispatch = old.auto), add = TRUE)
  on.exit(close_mpi_slaves(), add = TRUE)
  options(npRmpi.autodispatch = TRUE)

  set.seed(20260308)
  n <- 70
  x <- runif(n)
  z <- runif(n)
  y <- (0.4 + x) * sin(2 * pi * z) + rnorm(n, sd = 0.04)

  tx.sc <- data.frame(x = x)
  tz.sc <- data.frame(z = z)
  ex.sc <- data.frame(x = seq(0.1, 0.9, length.out = 18))
  ez.sc <- data.frame(z = seq(0.1, 0.9, length.out = 18))

  for (regtype in c("lc", "ll", "lp")) {
    for (bwtype in c("fixed", "generalized_nn", "adaptive_nn")) {
      sc.args <- list(
        xdat = tx.sc,
        zdat = tz.sc,
        ydat = y,
        regtype = regtype,
        bwtype = bwtype,
        bandwidth.compute = FALSE,
        bws = if (identical(bwtype, "fixed")) 0.18 else 9
      )
      if (identical(regtype, "lp")) {
        sc.args$degree <- 1L
        sc.args$basis <- "glp"
        sc.args$bernstein.basis <- FALSE
      }

      sc.bw <- do.call(npscoefbw, sc.args)
      sc.fit <- npscoef(
        bws = sc.bw,
        txdat = tx.sc,
        tzdat = tz.sc,
        tydat = y,
        exdat = ex.sc,
        ezdat = ez.sc,
        iterate = FALSE,
        errors = FALSE
      )
      sc.apply <- npscoefhat(
        bws = sc.bw,
        txdat = tx.sc,
        tzdat = tz.sc,
        exdat = ex.sc,
        ezdat = ez.sc,
        y = y,
        output = "apply",
        iterate = FALSE
      )
      sc.H <- npscoefhat(
        bws = sc.bw,
        txdat = tx.sc,
        tzdat = tz.sc,
        exdat = ex.sc,
        ezdat = ez.sc,
        output = "matrix",
        iterate = FALSE
      )
      expect_equal(
        as.vector(sc.apply),
        as.vector(sc.H %*% y),
        tolerance = 1e-10,
        info = paste("npscoefhat helper parity", regtype, bwtype)
      )
      if (identical(bwtype, "fixed")) {
        expect_equal(
          as.vector(sc.apply),
          as.vector(sc.fit$mean),
          tolerance = 1e-8,
          info = paste("npscoefhat core parity", regtype, bwtype)
        )
      }
    }
  }
})

test_that("npscoefhat selected adaptive-nn owner preserves integer support and helper parity", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  old.auto <- getOption("npRmpi.autodispatch", FALSE)
  on.exit(options(npRmpi.autodispatch = old.auto), add = TRUE)
  on.exit(close_mpi_slaves(), add = TRUE)
  options(npRmpi.autodispatch = TRUE)
  round_half_to_even <- getFromNamespace(".np_round_half_to_even", "npRmpi")

  set.seed(20260308)
  n <- 50L
  x <- runif(n)
  z <- runif(n)
  y <- (0.4 + x) * sin(2 * pi * z) + rnorm(n, sd = 0.04)

  tx <- data.frame(x = x)
  tz <- data.frame(z = z)
  ex <- data.frame(x = seq(0.1, 0.9, length.out = 12L))
  ez <- data.frame(z = seq(0.1, 0.9, length.out = 12L))

  tol <- sqrt(.Machine$double.eps)
  upper <- n - 1L

  for (regtype in c("lc", "ll", "lp")) {
    bw.args <- list(
      xdat = tx,
      zdat = tz,
      ydat = y,
      regtype = regtype,
      bwtype = "adaptive_nn",
      nmulti = 1L
    )
    if (identical(regtype, "lp")) {
      bw.args$degree <- 1L
      bw.args$basis <- "glp"
      bw.args$bernstein.basis <- FALSE
    }

    bw <- do.call(npscoefbw, bw.args)
    expect_true(all(abs(bw$bw - round_half_to_even(bw$bw)) <= tol), info = regtype)
    expect_true(all(bw$bw >= 1 & bw$bw <= upper), info = regtype)

    hat.apply <- npscoefhat(
      bws = bw,
      txdat = tx,
      tzdat = tz,
      exdat = ex,
      ezdat = ez,
      y = y,
      output = "apply",
      iterate = FALSE
    )
    hat.matrix <- npscoefhat(
      bws = bw,
      txdat = tx,
      tzdat = tz,
      exdat = ex,
      ezdat = ez,
      output = "matrix",
      iterate = FALSE
    )

    expect_equal(as.vector(hat.apply), as.vector(hat.matrix %*% y), tolerance = 1e-10, info = regtype)
  }
})


test_that("semihat validates class and scalar controls", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  old.auto <- getOption("npRmpi.autodispatch", FALSE)
  on.exit(options(npRmpi.autodispatch = old.auto), add = TRUE)
  on.exit(close_mpi_slaves(), add = TRUE)
  options(npRmpi.autodispatch = TRUE)

  set.seed(27182)
  n <- 40
  x <- runif(n)
  y <- rnorm(n)
  z <- runif(n)

  rbw <- npregbw(y ~ x, bws = 0.2, bandwidth.compute = FALSE)
  expect_error(npindexhat(bws = rbw, txdat = data.frame(x = x)), "sibandwidth")
  expect_error(npplreghat(bws = rbw, txdat = data.frame(x = x), tzdat = data.frame(z = z)), "plbandwidth")
  expect_error(npscoefhat(bws = rbw, txdat = data.frame(x = x), tzdat = data.frame(z = z)), "scbandwidth")

  sibw <- npindexbw(
    xdat = data.frame(x = x, x2 = x^2),
    ydat = y,
    bws = c(1, 1, 0.2),
    bandwidth.compute = FALSE,
    regtype = "ll"
  )
  expect_error(
    npindexhat(bws = sibw, txdat = data.frame(x = x, x2 = x^2), s = 1L, fd.step = 0),
    "argument 'fd.step' must be a positive finite scalar"
  )

  scbw <- npscoefbw(
    xdat = x,
    zdat = z,
    ydat = y,
    bws = 0.2,
    bandwidth.compute = FALSE
  )
  h0 <- npscoefhat(
    bws = scbw,
    txdat = data.frame(x = x),
    tzdat = data.frame(z = z),
    ridge = 0
  )
  expect_true(is.matrix(h0))
})

test_that("plot bootstrap supports wild for sc/pl/si bandwidth objects", {
  skip_on_cran()
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  old.auto <- getOption("npRmpi.autodispatch", FALSE)
  on.exit(options(npRmpi.autodispatch = old.auto), add = TRUE)
  on.exit(close_mpi_slaves(), add = TRUE)
  options(npRmpi.autodispatch = TRUE)

  old.chunk <- getOption("np.plot.wild.chunk.size")
  on.exit(options(np.plot.wild.chunk.size = old.chunk), add = TRUE)
  options(np.plot.wild.chunk.size = 5L)

  set.seed(20260223)
  n <- 80
  x <- runif(n)
  z <- runif(n)
  y <- sin(2 * pi * z) + 1.5 * x + rnorm(n, sd = 0.05)

  scbw <- npscoefbw(
    xdat = x,
    zdat = z,
    ydat = y,
    bws = 0.2,
    bandwidth.compute = FALSE
  )
  sc.out <- suppressWarnings(
    plot(
      scbw,
      xdat = data.frame(x = x),
      ydat = y,
      zdat = data.frame(z = z),
      perspective = FALSE,
      plot.behavior = "data",
      plot.errors.method = "bootstrap",
      plot.errors.boot.method = "wild",
      plot.errors.boot.num = 19
    )
  )
  expect_type(sc.out, "list")
  expect_true(length(sc.out) > 0)

  sc.out.rad <- suppressWarnings(
    plot(
      scbw,
      xdat = data.frame(x = x),
      ydat = y,
      zdat = data.frame(z = z),
      perspective = FALSE,
      plot.behavior = "data",
      plot.errors.method = "bootstrap",
      plot.errors.boot.method = "wild",
      plot.errors.boot.wild = "rademacher",
      plot.errors.type = "pointwise",
      plot.errors.boot.num = 19
    )
  )
  expect_type(sc.out.rad, "list")
  expect_true(length(sc.out.rad) > 0)
  expect_false(isTRUE(all.equal(
    sc.out[[1]]$merr,
    sc.out.rad[[1]]$merr,
    tolerance = 0,
    check.attributes = FALSE
  )))

  pl.cfgs <- list(
    list(regtype = "lc", basis = NULL, label = "lc"),
    list(regtype = "ll", basis = NULL, label = "ll"),
    list(regtype = "lp", basis = "glp", label = "lp-glp"),
    list(regtype = "lp", basis = "additive", label = "lp-additive"),
    list(regtype = "lp", basis = "tensor", label = "lp-tensor")
  )
  for (cfg in pl.cfgs) {
    pl.args <- list(
      xdat = x,
      zdat = z,
      ydat = y,
      bws = matrix(c(0.2, 0.2), nrow = 2, ncol = 1),
      bandwidth.compute = FALSE,
      regtype = cfg$regtype
    )
    if (!is.null(cfg$basis)) {
      pl.args$basis <- cfg$basis
      pl.args$degree <- 2L
    }
    plbw <- do.call(npplregbw, pl.args)
    pl.out <- suppressWarnings(
      plot(
        plbw,
        xdat = data.frame(x = x),
        ydat = y,
        zdat = data.frame(z = z),
        perspective = FALSE,
        plot.behavior = "data",
        plot.errors.method = "bootstrap",
        plot.errors.boot.method = "wild",
        plot.errors.boot.num = 9
      )
    )
    expect_true(is.list(pl.out), info = cfg$label)
    expect_true(length(pl.out) > 0, info = cfg$label)
  }

  sibw <- npindexbw(
    xdat = data.frame(x1 = x, x2 = z),
    ydat = y,
    method = "ichimura",
    nmulti = 1
  )
  si.out <- suppressWarnings(
    plot(
      sibw,
      xdat = data.frame(x1 = x, x2 = z),
      ydat = y,
      plot.behavior = "data",
      plot.errors.method = "bootstrap",
      plot.errors.boot.method = "wild",
      plot.errors.boot.num = 19
    )
  )
  expect_type(si.out, "list")
  expect_true(length(si.out) > 0)
  expect_true(is.matrix(si.out[[1]]$merr))
  expect_equal(ncol(si.out[[1]]$merr), 2L)
  expect_false(all(is.na(si.out[[1]]$merr)))
})
