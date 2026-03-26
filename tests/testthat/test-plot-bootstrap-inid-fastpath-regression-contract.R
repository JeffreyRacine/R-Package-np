test_that("npindex inid fast path matches explicit resample refits", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  old.auto <- getOption("npRmpi.autodispatch", FALSE)
  on.exit(options(npRmpi.autodispatch = old.auto), add = TRUE)
  on.exit(close_mpi_slaves(), add = TRUE)
  options(npRmpi.autodispatch = TRUE)

  set.seed(3231)
  n <- 40
  x1 <- runif(n)
  x2 <- runif(n)
  y <- sin(x1 + x2) + rnorm(n, sd = 0.1)
  tx <- data.frame(x1 = x1, x2 = x2)
  bw <- npindexbw(xdat = tx, ydat = y, method = "ichimura", nmulti = 1)
  B <- 11L
  counts <- rmultinom(n = B, size = n, prob = rep.int(1 / n, n))

  fast.fun <- getFromNamespace(".np_inid_lc_boot_from_hat", "npRmpi")
  H <- npindexhat(bws = bw, txdat = tx, exdat = tx, output = "matrix", s = 0L)
  fast.out <- fast.fun(H = H, ydat = y, B = B, counts = counts)

  explicit.t <- matrix(NA_real_, nrow = B, ncol = n)
  for (b in seq_len(B)) {
    idx <- rep.int(seq_len(n), counts[, b])
    explicit.t[b, ] <- npindex(
      txdat = tx[idx, , drop = FALSE],
      tydat = y[idx],
      exdat = tx,
      bws = bw,
      gradients = FALSE
    )$mean
  }

  expect_equal(fast.out$t, explicit.t, tolerance = 1e-10)
  expect_equal(fast.out$t0, as.vector(H %*% y), tolerance = 1e-12)
})

test_that("npindex ll/lp inid fast path matches explicit resample refits", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  old.auto <- getOption("npRmpi.autodispatch", FALSE)
  on.exit(options(npRmpi.autodispatch = old.auto), add = TRUE)
  on.exit(close_mpi_slaves(), add = TRUE)
  options(npRmpi.autodispatch = FALSE)

  set.seed(32315)
  n <- 45
  x1 <- runif(n)
  x2 <- runif(n)
  y <- sin(x1 + x2) + rnorm(n, sd = 0.1)
  tx <- data.frame(x1 = x1, x2 = x2)
  B <- 9L
  counts <- rmultinom(n = B, size = n, prob = rep.int(1 / n, n))

  fast.fun <- getFromNamespace(".np_inid_boot_from_regression", "npRmpi")
  rbw.fun <- getFromNamespace(".np_indexhat_rbw", "npRmpi")

  cfgs <- list(
    list(regtype = "ll", basis = NULL, label = "ll"),
    list(regtype = "lp", basis = "additive", label = "lp-additive"),
    list(regtype = "lp", basis = "tensor", label = "lp-tensor"),
    list(regtype = "lp", basis = "glp", label = "lp-glp")
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

    tx.index <- data.frame(index = as.vector(as.matrix(tx) %*% bw$beta))
    rbw <- rbw.fun(bws = bw, idx.train = tx.index)

    fast.out <- fast.fun(
      xdat = tx.index,
      exdat = tx.index,
      bws = rbw,
      ydat = y,
      B = B,
      counts = counts
    )

    explicit.t <- matrix(NA_real_, nrow = B, ncol = n)
    for (b in seq_len(B)) {
      idx <- rep.int(seq_len(n), counts[, b])
      explicit.t[b, ] <- npindex(
        txdat = tx[idx, , drop = FALSE],
        tydat = y[idx],
        exdat = tx,
        bws = bw,
        gradients = FALSE
      )$mean
    }

    expect_equal(fast.out$t, explicit.t, tolerance = 1e-6, info = cfg$label)
    expect_equal(
      fast.out$t0,
      npindex(txdat = tx, tydat = y, exdat = tx, bws = bw, gradients = FALSE)$mean,
      tolerance = 1e-6,
      info = cfg$label
    )
  }
})

test_that("npreg inid ll/lp fast path matches explicit resample refits", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  old.auto <- getOption("npRmpi.autodispatch", FALSE)
  on.exit(options(npRmpi.autodispatch = old.auto), add = TRUE)
  on.exit(close_mpi_slaves(), add = TRUE)
  options(npRmpi.autodispatch = FALSE)

  set.seed(3211)
  n <- 45
  x1 <- runif(n)
  x2 <- runif(n)
  y <- sin(2 * pi * x1) + 0.6 * x2 + rnorm(n, sd = 0.08)
  tx <- data.frame(x1 = x1, x2 = x2)
  ex <- tx[seq_len(15), , drop = FALSE]
  B <- 9L
  counts <- rmultinom(n = B, size = n, prob = rep.int(1 / n, n))

  fast.fun <- getFromNamespace(".np_inid_boot_from_regression", "npRmpi")
  cfgs <- list(
    list(regtype = "ll", basis = NULL, degree = NULL, label = "ll"),
    list(regtype = "lp", basis = "additive", degree = c(2L, 2L), label = "lp-additive"),
    list(regtype = "lp", basis = "tensor", degree = c(2L, 2L), label = "lp-tensor"),
    list(regtype = "lp", basis = "glp", degree = c(2L, 2L), label = "lp-glp")
  )

  for (cfg in cfgs) {
    bw.args <- list(
      xdat = tx,
      ydat = y,
      regtype = cfg$regtype,
      bws = c(0.3, 0.3),
      bandwidth.compute = FALSE
    )
    if (!is.null(cfg$basis)) {
      bw.args$basis <- cfg$basis
      bw.args$degree <- cfg$degree
    }
    bw <- do.call(npregbw, bw.args)

    fast.out <- fast.fun(
      xdat = tx,
      exdat = ex,
      bws = bw,
      ydat = y,
      B = B,
      counts = counts
    )

    explicit.t <- matrix(NA_real_, nrow = B, ncol = nrow(ex))
    for (b in seq_len(B)) {
      idx <- rep.int(seq_len(n), counts[, b])
      explicit.t[b, ] <- npreg(
        txdat = tx[idx, , drop = FALSE],
        tydat = y[idx],
        exdat = ex,
        bws = bw,
        gradients = FALSE,
        warn.glp.gradient = FALSE
      )$mean
    }

    expect_equal(fast.out$t, explicit.t, tolerance = 1e-6, info = cfg$label)
    expect_equal(
      fast.out$t0,
      npreg(txdat = tx, tydat = y, exdat = ex, bws = bw, gradients = FALSE, warn.glp.gradient = FALSE)$mean,
      tolerance = 1e-6,
      info = cfg$label
    )
  }
})

test_that("npindex bounded bootstrap plot-data supports bw and fit objects", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  old.auto <- getOption("npRmpi.autodispatch", FALSE)
  on.exit(options(npRmpi.autodispatch = old.auto), add = TRUE)
  on.exit(close_mpi_slaves(), add = TRUE)
  options(npRmpi.autodispatch = TRUE)

  set.seed(32322)
  n <- 36
  x1 <- runif(n)
  x2 <- runif(n)
  y <- sin(x1 + x2) + rnorm(n, sd = 0.08)
  tx <- data.frame(x1 = x1, x2 = x2)

  bw <- npindexbw(
    xdat = tx,
    ydat = y,
    method = "ichimura",
    ckerbound = "range",
    bws = c(1, 0.25, 0.25),
    bandwidth.compute = FALSE
  )
  fit <- npindex(bws = bw, txdat = tx, tydat = y)

  run_plot <- function(obj, boot.method) {
    suppressWarnings(
      plot(
        obj,
        xdat = tx,
        ydat = y,
        plot.behavior = "data",
        perspective = FALSE,
        plot.errors.method = "bootstrap",
        plot.errors.boot.method = boot.method,
        plot.errors.boot.num = 5L,
        plot.errors.type = "pointwise"
      )
    )
  }

  for (boot.method in c("inid", "geom", "wild")) {
    out.bw <- run_plot(bw, boot.method)
    out.fit <- run_plot(fit, boot.method)

    expect_type(out.bw, "list")
    expect_type(out.fit, "list")
    expect_true(length(out.bw) > 0L, info = paste("bw", boot.method))
    expect_true(length(out.fit) > 0L, info = paste("fit", boot.method))
  }
})

test_that("npplreg inid fast path matches explicit resample refits", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  old.auto <- getOption("npRmpi.autodispatch", FALSE)
  on.exit(options(npRmpi.autodispatch = old.auto), add = TRUE)
  on.exit(close_mpi_slaves(), add = TRUE)
  options(npRmpi.autodispatch = TRUE)

  set.seed(32316)
  n <- 40
  x1 <- runif(n)
  x2 <- runif(n)
  z1 <- runif(n)
  z2 <- runif(n)
  tx <- data.frame(x1 = x1, x2 = x2)
  tz <- data.frame(z1 = z1, z2 = z2)
  y <- sin(2 * pi * z1) + 0.5 * x1 - 0.2 * x2 + rnorm(n, sd = 0.08)
  B <- 9L
  counts <- rmultinom(n = B, size = n, prob = rep.int(1 / n, n))

  bw <- npplregbw(xdat = tx, ydat = y, zdat = tz, regtype = "lc", nmulti = 1)
  fast.fun <- getFromNamespace(".np_inid_boot_from_plreg", "npRmpi")
  fast.out <- fast.fun(
    txdat = tx,
    ydat = y,
    tzdat = tz,
    exdat = tx,
    ezdat = tz,
    bws = bw,
    B = B,
    counts = counts
  )

  explicit.t <- matrix(NA_real_, nrow = B, ncol = n)
  for (b in seq_len(B)) {
    idx <- rep.int(seq_len(n), counts[, b])
    explicit.t[b, ] <- npplreg(
      bws = bw,
      txdat = tx[idx, , drop = FALSE],
      tydat = y[idx],
      tzdat = tz[idx, , drop = FALSE],
      exdat = tx,
      ezdat = tz
    )$mean
  }

  fit0 <- npplreg(
    bws = bw,
    txdat = tx,
    tydat = y,
    tzdat = tz,
    exdat = tx,
    ezdat = tz
  )$mean

  expect_equal(fast.out$t, explicit.t, tolerance = 1e-6)
  expect_equal(as.vector(fast.out$t0), as.vector(fit0), tolerance = 1e-7)
})

test_that("npreg inid fast path supports continuous-slice gradients", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  old.auto <- getOption("npRmpi.autodispatch", FALSE)
  on.exit(options(npRmpi.autodispatch = old.auto), add = TRUE)
  on.exit(close_mpi_slaves(), add = TRUE)
  options(npRmpi.autodispatch = FALSE)

  set.seed(3212)
  n <- 45
  x1 <- runif(n)
  x2 <- runif(n)
  y <- sin(2 * pi * x1) + 0.4 * x2 + rnorm(n, sd = 0.08)
  tx <- data.frame(x1 = x1, x2 = x2)
  ex <- tx[seq_len(15), , drop = FALSE]
  B <- 9L
  counts <- rmultinom(n = B, size = n, prob = rep.int(1 / n, n))

  bw <- npregbw(
    xdat = tx,
    ydat = y,
    regtype = "ll",
    bws = c(0.3, 0.3),
    bandwidth.compute = FALSE
  )

  fast.fun <- getFromNamespace(".np_inid_boot_from_regression", "npRmpi")
  fast.out <- fast.fun(
    xdat = tx,
    exdat = ex,
    bws = bw,
    ydat = y,
    B = B,
    counts = counts,
    gradients = TRUE,
    gradient.order = 1L,
    slice.index = 1L
  )

  explicit.t <- matrix(NA_real_, nrow = B, ncol = nrow(ex))
  for (b in seq_len(B)) {
    idx <- rep.int(seq_len(n), counts[, b])
    explicit.t[b, ] <- npreg(
      txdat = tx[idx, , drop = FALSE],
      tydat = y[idx],
      exdat = ex,
      bws = bw,
      gradients = TRUE,
      gradient.order = 1L,
      warn.glp.gradient = FALSE
    )$grad[, 1L]
  }

  fit0 <- npreg(
    txdat = tx,
    tydat = y,
    exdat = ex,
    bws = bw,
    gradients = TRUE,
    gradient.order = 1L,
    warn.glp.gradient = FALSE
  )$grad[, 1L]

  expect_equal(fast.out$t, explicit.t, tolerance = 1e-6)
  expect_equal(as.vector(fast.out$t0), as.vector(fit0), tolerance = 1e-6)
})

test_that("rbandwidth plot bootstrap supports gradients across methods", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  old.auto <- getOption("npRmpi.autodispatch", FALSE)
  on.exit(options(npRmpi.autodispatch = old.auto), add = TRUE)
  on.exit(close_mpi_slaves(), add = TRUE)
  options(npRmpi.autodispatch = TRUE)

  set.seed(322)
  n <- 60
  x <- runif(n)
  y <- sin(2 * pi * x) + rnorm(n, sd = 0.15)
  bw <- npregbw(y ~ x, regtype = "ll", nmulti = 1)

  methods <- c("inid", "fixed", "geom", "wild")
  for (m in methods) {
    set.seed(8320 + match(m, methods))
    out <- suppressWarnings(
      plot(
        bw,
        plot.behavior = "data",
        perspective = FALSE,
        gradients = TRUE,
        plot.errors.method = "bootstrap",
        plot.errors.boot.method = m,
        plot.errors.boot.num = 9
      )
    )
    expect_type(out, "list")
    expect_true(length(out) > 0)
  }
})
