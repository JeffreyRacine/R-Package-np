test_that("inid lc fast path matches explicit resample refits", {
  skip_if_not_installed("np")

  set.seed(321)
  n <- 40
  x <- runif(n)
  y <- cos(2 * pi * x) + rnorm(n, sd = 0.1)
  tx <- data.frame(x = x)
  ex <- data.frame(x = seq(min(x), max(x), length.out = 21))

  bw <- npregbw(
    y ~ x,
    regtype = "lc",
    bws = c(0.2),
    bandwidth.compute = FALSE
  )

  H <- npreghat(bws = bw, txdat = tx, exdat = ex, output = "matrix")
  B <- 13L
  counts <- rmultinom(n = B, size = n, prob = rep.int(1 / n, n))

  fast.fun <- getFromNamespace(".np_inid_lc_boot_from_hat", "np")
  fast.out <- fast.fun(H = H, ydat = y, B = B, counts = counts)

  explicit.t <- matrix(NA_real_, nrow = B, ncol = nrow(ex))
  for (b in seq_len(B)) {
    idx <- rep.int(seq_len(n), counts[, b])
    fit.b <- npreg(
      txdat = tx[idx, , drop = FALSE],
      tydat = y[idx],
      exdat = ex,
      bws = bw,
      gradients = FALSE,
      warn.glp.gradient = FALSE
    )
    explicit.t[b, ] <- fit.b$mean
  }

  expect_equal(fast.out$t, explicit.t, tolerance = 1e-10)
  expect_equal(fast.out$t0, as.vector(H %*% y), tolerance = 1e-12)
})

test_that("inid ll/lp fast path matches explicit resample refits", {
  skip_if_not_installed("np")

  set.seed(3211)
  n <- 45
  x1 <- runif(n)
  x2 <- runif(n)
  y <- sin(2 * pi * x1) + 0.6 * x2 + rnorm(n, sd = 0.08)
  tx <- data.frame(x1 = x1, x2 = x2)
  ex <- tx[seq_len(15), , drop = FALSE]
  B <- 9L
  counts <- rmultinom(n = B, size = n, prob = rep.int(1 / n, n))

  fast.fun <- getFromNamespace(".np_inid_boot_from_regression", "np")
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
  skip_if_not_installed("np")

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
  fast.fun <- getFromNamespace(".np_inid_boot_from_plreg", "np")
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
  skip_if_not_installed("np")

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

  fast.fun <- getFromNamespace(".np_inid_boot_from_regression", "np")
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
