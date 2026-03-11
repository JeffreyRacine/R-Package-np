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

test_that("npscoef plot bootstrap inid supports ll/lp basis variants", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  old.auto <- getOption("npRmpi.autodispatch", FALSE)
  on.exit(options(npRmpi.autodispatch = old.auto), add = TRUE)
  on.exit(close_mpi_slaves(), add = TRUE)
  options(npRmpi.autodispatch = TRUE)

  set.seed(3233)
  n <- 70
  x <- runif(n)
  z <- runif(n)
  y <- (0.6 + x) * sin(2 * pi * z) + rnorm(n, sd = 0.08)
  tx <- data.frame(x = x)
  tz <- data.frame(z = z)

  cfgs <- list(
    list(regtype = "ll", basis = NULL, label = "ll"),
    list(regtype = "lp", basis = "glp", label = "lp-glp"),
    list(regtype = "lp", basis = "additive", label = "lp-additive"),
    list(regtype = "lp", basis = "tensor", label = "lp-tensor")
  )

  for (cfg in cfgs) {
    bw.args <- list(
      xdat = tx,
      zdat = tz,
      ydat = y,
      bws = 0.22,
      bandwidth.compute = FALSE,
      regtype = cfg$regtype
    )
    if (!is.null(cfg$basis)) {
      bw.args$basis <- cfg$basis
      bw.args$degree <- 2L
    }
    bw <- do.call(npscoefbw, bw.args)
    out <- suppressWarnings(
      plot(
        bw,
        xdat = tx,
        ydat = y,
        zdat = tz,
        perspective = FALSE,
        plot.behavior = "data",
        plot.errors.method = "bootstrap",
        plot.errors.boot.method = "inid",
        plot.errors.boot.num = 7
      )
    )
    expect_type(out, "list")
    expect_true(length(out) > 0, info = cfg$label)
  }
})

test_that("density/distribution plot bootstrap rejects wild selector", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  old.auto <- getOption("npRmpi.autodispatch", FALSE)
  on.exit(options(npRmpi.autodispatch = old.auto), add = TRUE)
  on.exit(close_mpi_slaves(), add = TRUE)
  options(npRmpi.autodispatch = TRUE)

  set.seed(328)
  n <- 50
  x <- rnorm(n)
  y <- rnorm(n)

  ubw <- npudensbw(dat = data.frame(x = x), bws = 0.5, bandwidth.compute = FALSE)
  expect_error(
    suppressWarnings(plot(
      ubw,
      plot.behavior = "data",
      perspective = FALSE,
      plot.errors.method = "bootstrap",
      plot.errors.boot.method = "wild",
      plot.errors.boot.num = 9
    )),
    "not supported for unconditional density/distribution estimators"
  )

  cbw <- npcdensbw(
    xdat = data.frame(x = x),
    ydat = data.frame(y = y),
    bws = c(0.6, 0.6),
    bandwidth.compute = FALSE
  )
  expect_error(
    suppressWarnings(plot(
      cbw,
      plot.behavior = "data",
      perspective = FALSE,
      plot.errors.method = "bootstrap",
      plot.errors.boot.method = "wild",
      plot.errors.boot.num = 9
    )),
    "not supported for conditional density/distribution estimators"
  )
})

test_that("nearest-neighbor plot helpers run for regression and density/distribution families", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  old.auto <- getOption("npRmpi.autodispatch", FALSE)
  on.exit(options(npRmpi.autodispatch = old.auto), add = TRUE)
  on.exit(close_mpi_slaves(), add = TRUE)
  options(npRmpi.autodispatch = TRUE)

  set.seed(329)
  n <- 45
  x <- data.frame(x = runif(n))
  y <- sin(2 * pi * x$x) + rnorm(n, sd = 0.08)
  yframe <- data.frame(y = y)

  run_plot <- function(bw, ...) {
    suppressWarnings(plot(
      bw,
      plot.behavior = "data",
      perspective = FALSE,
      plot.errors.method = "bootstrap",
      plot.errors.boot.method = "inid",
      plot.errors.boot.num = 5,
      ...
    ))
  }

  for (bt in c("generalized_nn", "adaptive_nn")) {
    rbw <- npregbw(xdat = x, ydat = y, regtype = "ll", nmulti = 1, bwtype = bt)
    expect_type(run_plot(rbw, xdat = x, ydat = y), "list")
    expect_type(run_plot(rbw, xdat = x, ydat = y, gradients = TRUE), "list")

    ubw <- npudensbw(dat = x, nmulti = 1, bwtype = bt)
    expect_type(run_plot(ubw), "list")

    dbw <- npudistbw(dat = x, nmulti = 1, bwtype = bt)
    expect_type(run_plot(dbw), "list")

    cbw <- npcdensbw(xdat = x, ydat = yframe, nmulti = 1, bwtype = bt)
    expect_type(run_plot(cbw, xdat = x, ydat = yframe, view = "fixed"), "list")

    cdbw <- npcdistbw(xdat = x, ydat = yframe, nmulti = 1, bwtype = bt)
    expect_type(run_plot(cdbw, xdat = x, ydat = yframe, view = "fixed"), "list")
  }
})

test_that("nearest-neighbor regression and ksum helpers match explicit resample refits", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  old.auto <- getOption("npRmpi.autodispatch", FALSE)
  on.exit(options(npRmpi.autodispatch = old.auto), add = TRUE)
  on.exit(close_mpi_slaves(), add = TRUE)
  options(npRmpi.autodispatch = FALSE)

  set.seed(32905)
  n <- 40
  x1 <- runif(n)
  x2 <- runif(n)
  y <- sin(2 * pi * x1) + 0.5 * x2 + rnorm(n, sd = 0.08)
  tx <- data.frame(x1 = x1, x2 = x2)
  ex <- tx[seq_len(12), , drop = FALSE]
  ty <- data.frame(y = runif(n))
  ey <- data.frame(y = seq(0.1, 0.9, length.out = 12))
  B <- 7L
  counts <- rmultinom(n = B, size = n, prob = rep.int(1 / n, n))

  reg.fast <- getFromNamespace(".np_inid_boot_from_regression", "npRmpi")
  fast.u <- getFromNamespace(".np_inid_boot_from_ksum_unconditional", "npRmpi")
  fast.c <- getFromNamespace(".np_inid_boot_from_ksum_conditional", "npRmpi")

  for (bt in c("generalized_nn", "adaptive_nn")) {
    bw.val <- if (identical(bt, "adaptive_nn")) c(5, 5) else c(1, 1)
    rbw <- npregbw(
      xdat = tx,
      ydat = y,
      regtype = "ll",
      bws = bw.val,
      bwtype = bt,
      bandwidth.compute = FALSE
    )
    reg.out <- reg.fast(
      xdat = tx,
      exdat = ex,
      bws = rbw,
      ydat = y,
      B = B,
      counts = counts
    )
    reg.explicit <- matrix(NA_real_, nrow = B, ncol = nrow(ex))
    for (b in seq_len(B)) {
      idx <- rep.int(seq_len(n), counts[, b])
      reg.explicit[b, ] <- npreg(
        txdat = tx[idx, , drop = FALSE],
        tydat = y[idx],
        exdat = ex,
        bws = rbw,
        gradients = FALSE,
        warn.glp.gradient = FALSE
      )$mean
    }
    expect_equal(reg.out$t, reg.explicit, tolerance = 1e-10, info = paste(bt, "npreg"))
    expect_equal(
      reg.out$t0,
      npreg(txdat = tx, tydat = y, exdat = ex, bws = rbw, gradients = FALSE, warn.glp.gradient = FALSE)$mean,
      tolerance = 1e-10,
      info = paste(bt, "npreg-t0")
    )

    ubw <- npudensbw(dat = data.frame(x = x1), bws = bw.val[1L], bwtype = bt, bandwidth.compute = FALSE)
    uout <- fast.u(
      xdat = data.frame(x = x1),
      exdat = data.frame(x = ex$x1),
      bws = ubw,
      B = B,
      operator = "normal",
      counts = counts
    )
    uexplicit <- matrix(NA_real_, nrow = B, ncol = nrow(ex))
    for (b in seq_len(B)) {
      idx <- rep.int(seq_len(n), counts[, b])
      uexplicit[b, ] <- npudens(
        tdat = data.frame(x = x1)[idx, , drop = FALSE],
        edat = data.frame(x = ex$x1),
        bws = ubw
      )$dens
    }
    expect_equal(uout$t, uexplicit, tolerance = 1e-10, info = paste(bt, "npudens"))

    cbw <- npcdensbw(
      xdat = data.frame(x = x1),
      ydat = ty,
      bws = rep.int(bw.val[1L], 2L),
      bwtype = bt,
      bandwidth.compute = FALSE
    )
    cout <- fast.c(
      xdat = data.frame(x = x1),
      ydat = ty,
      exdat = data.frame(x = ex$x1),
      eydat = ey,
      bws = cbw,
      B = B,
      cdf = FALSE,
      counts = counts
    )
    cexplicit <- matrix(NA_real_, nrow = B, ncol = nrow(ex))
    for (b in seq_len(B)) {
      idx <- rep.int(seq_len(n), counts[, b])
      cexplicit[b, ] <- npcdens(
        txdat = data.frame(x = x1)[idx, , drop = FALSE],
        tydat = ty[idx, , drop = FALSE],
        exdat = data.frame(x = ex$x1),
        eydat = ey,
        bws = cbw
      )$condens
    }
    expect_equal(cout$t, cexplicit, tolerance = 1e-10, info = paste(bt, "npcdens"))
  }
})
