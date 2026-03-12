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

test_that("rbandwidth plot bootstrap supports gradients across methods", {
  skip_if_not_installed("np")

  set.seed(322)
  n <- 60
  x <- runif(n)
  y <- sin(2 * pi * x) + rnorm(n, sd = 0.15)
  bw <- npregbw(y ~ x, regtype = "ll", nmulti = 1)

  run_plot <- function(method) {
    suppressWarnings(
      plot(
        bw,
        plot.behavior = "data",
        perspective = FALSE,
        gradients = TRUE,
        plot.errors.method = "bootstrap",
        plot.errors.boot.method = method,
        plot.errors.boot.num = 9
      )
    )
  }

  methods <- c("inid", "fixed", "geom", "wild")
  for (m in methods) {
    set.seed(9320 + match(m, methods))
    out <- run_plot(method = m)
    expect_type(out, "list")
    expect_true(length(out) > 0)
  }
})

test_that("plot bootstrap accepts wild selector", {
  skip_if_not_installed("np")

  set.seed(323)
  n <- 80
  x <- runif(n)
  y <- sin(2 * pi * x) + rnorm(n, sd = 0.1)
  bw <- npregbw(y ~ x, bws = 0.2, bandwidth.compute = FALSE)

  out <- suppressWarnings(
    plot(
      bw,
      xdat = data.frame(x = x),
      ydat = y,
      plot.behavior = "data",
      perspective = FALSE,
      plot.errors.method = "bootstrap",
      plot.errors.boot.method = "wild",
      plot.errors.boot.num = 9
    )
  )

  expect_type(out, "list")
  expect_true(length(out) > 0)
})

test_that("scbandwidth bootstrap non-coef path avoids npscoef refits", {
  skip_if_not_installed("np")

  set.seed(3240)
  n <- 65
  x <- runif(n)
  z <- runif(n)
  y <- sin(2 * pi * z) + x * (1 + z) + rnorm(n, sd = 0.1)
  xdat <- data.frame(x = x)
  zdat <- data.frame(z = z)
  bw <- npscoefbw(xdat = xdat, ydat = y, zdat = zdat, regtype = "lc", nmulti = 1)

  np.ns <- asNamespace("np")
  ctr <- new.env(parent = emptyenv())
  ctr$n <- 0L
  trace(
    what = "npscoef",
    where = np.ns,
    tracer = bquote(.(ctr)$n <- .(ctr)$n + 1L),
    print = FALSE
  )
  on.exit(untrace("npscoef", where = np.ns), add = TRUE)

  out <- suppressWarnings(
    plot(
      bw,
      xdat = xdat,
      ydat = y,
      zdat = zdat,
      coef = FALSE,
      perspective = FALSE,
      plot.behavior = "data",
      plot.errors.method = "bootstrap",
      plot.errors.boot.method = "inid",
      plot.errors.boot.num = 7
    )
  )

  expect_type(out, "list")
  expect_true(length(out) > 0)
  expect_identical(ctr$n, 0L)
})

test_that("npindex inid fast path matches explicit resample refits", {
  skip_if_not_installed("np")

  set.seed(3231)
  n <- 40
  x1 <- runif(n)
  x2 <- runif(n)
  y <- sin(x1 + x2) + rnorm(n, sd = 0.1)
  tx <- data.frame(x1 = x1, x2 = x2)
  bw <- npindexbw(xdat = tx, ydat = y, method = "ichimura", nmulti = 1)
  B <- 11L
  counts <- rmultinom(n = B, size = n, prob = rep.int(1 / n, n))

  fast.fun <- getFromNamespace(".np_inid_lc_boot_from_hat", "np")
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
  skip_if_not_installed("np")

  set.seed(32315)
  n <- 45
  x1 <- runif(n)
  x2 <- runif(n)
  y <- sin(x1 + x2) + rnorm(n, sd = 0.1)
  tx <- data.frame(x1 = x1, x2 = x2)
  B <- 9L
  counts <- rmultinom(n = B, size = n, prob = rep.int(1 / n, n))

  fast.fun <- getFromNamespace(".np_inid_boot_from_regression", "np")
  rbw.fun <- getFromNamespace(".np_indexhat_rbw", "np")

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

test_that("npindex plot bootstrap inid supports ll/lp basis variants", {
  skip_if_not_installed("np")

  set.seed(3232)
  n <- 70
  x1 <- runif(n)
  x2 <- runif(n)
  y <- sin(x1 + x2) + rnorm(n, sd = 0.1)
  tx <- data.frame(x1 = x1, x2 = x2)

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
    out <- suppressWarnings(
      plot(
        bw,
        xdat = tx,
        ydat = y,
        plot.behavior = "data",
        perspective = FALSE,
        plot.errors.method = "bootstrap",
        plot.errors.boot.method = "inid",
        plot.errors.boot.num = 7
      )
    )
    expect_type(out, "list")
    expect_true(length(out) > 0, info = cfg$label)
  }
})

test_that("npindex plot bootstrap inid fails fast for unsupported nonfixed gradients", {
  skip_if_not_installed("np")

  set.seed(32321)
  n <- 60
  x1 <- runif(n)
  x2 <- runif(n)
  y <- sin(x1 + x2) + rnorm(n, sd = 0.1)
  tx <- data.frame(x1 = x1, x2 = x2)

  for (bt in c("generalized_nn", "adaptive_nn")) {
    bw <- npindexbw(
      xdat = tx,
      ydat = y,
      bws = c(1, 1, 5L),
      bandwidth.compute = FALSE,
      bwtype = bt
    )

    expect_error(
      suppressWarnings(
        plot(
          bw,
          xdat = tx,
          ydat = y,
          plot.behavior = "data",
          perspective = FALSE,
          gradients = TRUE,
          plot.errors.method = "bootstrap",
          plot.errors.boot.method = "inid",
          plot.errors.boot.num = 7
        )
      ),
      "gradients=FALSE",
      info = bt
    )
  }
})

test_that("npscoef plot bootstrap inid supports ll/lp basis variants", {
  skip_if_not_installed("np")

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

test_that("inid lc chunked generation matches explicit chunked draws", {
  skip_if_not_installed("np")

  set.seed(324)
  n <- 30
  x <- runif(n)
  y <- cos(2 * pi * x) + rnorm(n, sd = 0.1)
  tx <- data.frame(x = x)
  ex <- data.frame(x = seq(min(x), max(x), length.out = 17))
  bw <- npregbw(y ~ x, regtype = "lc", bws = c(0.2), bandwidth.compute = FALSE)
  H <- npreghat(bws = bw, txdat = tx, exdat = ex, output = "matrix")
  B <- 11L

  fast.fun <- getFromNamespace(".np_inid_lc_boot_from_hat", "np")
  old.chunk <- getOption("np.plot.inid.chunk.size")
  on.exit(options(np.plot.inid.chunk.size = old.chunk), add = TRUE)
  options(np.plot.inid.chunk.size = 3L)

  set.seed(7001)
  auto <- fast.fun(H = H, ydat = y, B = B, counts = NULL)

  set.seed(7001)
  counts <- matrix(NA_integer_, nrow = n, ncol = B)
  start <- 1L
  while (start <= B) {
    stopi <- min(B, start + 3L - 1L)
    bsz <- stopi - start + 1L
    counts[, start:stopi] <- rmultinom(n = bsz, size = n, prob = rep.int(1 / n, n))
    start <- stopi + 1L
  }
  explicit <- fast.fun(H = H, ydat = y, B = B, counts = counts)

  expect_equal(auto$t, explicit$t, tolerance = 1e-12)
  expect_equal(auto$t0, explicit$t0, tolerance = 1e-12)
})

test_that("inid ksum fast path matches explicit resample refits for npudens/npudist", {
  skip_if_not_installed("np")

  set.seed(325)
  n <- 70
  x <- rnorm(n)
  tx <- data.frame(x = x)
  ex <- data.frame(x = seq(min(x), max(x), length.out = 19))
  B <- 9L
  counts <- rmultinom(n = B, size = n, prob = rep.int(1 / n, n))

  dens.bw <- npudensbw(dat = tx, bws = 0.35, bandwidth.compute = FALSE)
  dist.bw <- npudistbw(dat = tx, bws = 0.35, bandwidth.compute = FALSE)
  fast.fun <- getFromNamespace(".np_inid_boot_from_ksum_unconditional", "np")

  dens.fast <- fast.fun(
    xdat = tx,
    exdat = ex,
    bws = dens.bw,
    B = B,
    operator = "normal",
    counts = counts
  )
  dist.fast <- fast.fun(
    xdat = tx,
    exdat = ex,
    bws = dist.bw,
    B = B,
    operator = "integral",
    counts = counts
  )

  dens.explicit <- matrix(NA_real_, nrow = B, ncol = nrow(ex))
  dist.explicit <- matrix(NA_real_, nrow = B, ncol = nrow(ex))
  for (b in seq_len(B)) {
    idx <- rep.int(seq_len(n), counts[, b])
    dens.explicit[b, ] <- npudens(tdat = tx[idx, , drop = FALSE], edat = ex, bws = dens.bw)$dens
    dist.explicit[b, ] <- npudist(tdat = tx[idx, , drop = FALSE], edat = ex, bws = dist.bw)$dist
  }

  expect_equal(dens.fast$t, dens.explicit, tolerance = 1e-10)
  expect_equal(dist.fast$t, dist.explicit, tolerance = 1e-10)
  expect_equal(dens.fast$t0, npudens(tdat = tx, edat = ex, bws = dens.bw)$dens, tolerance = 1e-12)
  expect_equal(dist.fast$t0, npudist(tdat = tx, edat = ex, bws = dist.bw)$dist, tolerance = 1e-12)
})

test_that("inid ksum fast path matches explicit resample refits for npcdens/npcdist", {
  skip_if_not_installed("np")

  set.seed(326)
  n <- 80
  x <- rnorm(n)
  y <- rnorm(n)
  tx <- data.frame(x = x)
  ty <- data.frame(y = y)
  ex <- data.frame(x = seq(min(x), max(x), length.out = 17))
  ey <- data.frame(y = seq(min(y), max(y), length.out = 17))
  B <- 7L
  counts <- rmultinom(n = B, size = n, prob = rep.int(1 / n, n))

  dens.bw <- npcdensbw(xdat = tx, ydat = ty, bws = c(0.45, 0.45), bandwidth.compute = FALSE)
  dist.bw <- npcdistbw(xdat = tx, ydat = ty, bws = c(0.45, 0.45), bandwidth.compute = FALSE)
  fast.fun <- getFromNamespace(".np_inid_boot_from_ksum_conditional", "np")

  dens.fast <- fast.fun(
    xdat = tx,
    ydat = ty,
    exdat = ex,
    eydat = ey,
    bws = dens.bw,
    B = B,
    cdf = FALSE,
    counts = counts
  )
  dist.fast <- fast.fun(
    xdat = tx,
    ydat = ty,
    exdat = ex,
    eydat = ey,
    bws = dist.bw,
    B = B,
    cdf = TRUE,
    counts = counts
  )

  expect_true(is.list(dens.fast))
  expect_true(is.list(dist.fast))

  dens.explicit <- matrix(NA_real_, nrow = B, ncol = nrow(ex))
  dist.explicit <- matrix(NA_real_, nrow = B, ncol = nrow(ex))
  for (b in seq_len(B)) {
    idx <- rep.int(seq_len(n), counts[, b])
    dens.explicit[b, ] <- npcdens(
      txdat = tx[idx, , drop = FALSE],
      tydat = ty[idx, , drop = FALSE],
      exdat = ex,
      eydat = ey,
      bws = dens.bw
    )$condens
    dist.explicit[b, ] <- npcdist(
      txdat = tx[idx, , drop = FALSE],
      tydat = ty[idx, , drop = FALSE],
      exdat = ex,
      eydat = ey,
      bws = dist.bw
    )$condist
  }

  expect_equal(dens.fast$t, dens.explicit, tolerance = 1e-10)
  expect_equal(dist.fast$t, dist.explicit, tolerance = 1e-10)
  expect_equal(dens.fast$t0, npcdens(txdat = tx, tydat = ty, exdat = ex, eydat = ey, bws = dens.bw)$condens, tolerance = 1e-12)
  expect_equal(dist.fast$t0, npcdist(txdat = tx, tydat = ty, exdat = ex, eydat = ey, bws = dist.bw)$condist, tolerance = 1e-12)
})

test_that("ksum fast paths honor non-default kernel/bound options for density/distribution families", {
  skip_if_not_installed("np")

  set.seed(3261)
  n <- 60
  x <- runif(n)
  y <- runif(n)
  tx <- data.frame(x = x)
  ty <- data.frame(y = y)
  ex <- data.frame(x = seq(0.02, 0.98, length.out = 13))
  ey <- data.frame(y = seq(0.03, 0.97, length.out = 13))
  B <- 6L
  counts <- rmultinom(n = B, size = n, prob = rep.int(1 / n, n))

  u.dens.bw <- npudensbw(
    dat = tx,
    bws = 0.20,
    bandwidth.compute = FALSE,
    bwtype = "fixed",
    ckertype = "epanechnikov",
    ckerbound = "fixed",
    ckerlb = 0.0,
    ckerub = 1.0
  )
  u.dist.bw <- npudistbw(
    dat = tx,
    bws = 0.20,
    bandwidth.compute = FALSE,
    bwtype = "fixed",
    ckertype = "epanechnikov",
    ckerbound = "fixed",
    ckerlb = 0.0,
    ckerub = 1.0
  )
  c.dens.bw <- npcdensbw(
    xdat = tx,
    ydat = ty,
    bws = c(0.25, 0.25),
    bandwidth.compute = FALSE,
    bwtype = "fixed",
    cxkertype = "epanechnikov",
    cykertype = "epanechnikov",
    cxkerbound = "fixed",
    cykerbound = "fixed",
    cxkerlb = 0.0,
    cykerlb = 0.0,
    cxkerub = 1.0,
    cykerub = 1.0
  )
  c.dist.bw <- npcdistbw(
    xdat = tx,
    ydat = ty,
    bws = c(0.25, 0.25),
    bandwidth.compute = FALSE,
    bwtype = "fixed",
    cxkertype = "epanechnikov",
    cykertype = "epanechnikov",
    cxkerbound = "fixed",
    cykerbound = "fixed",
    cxkerlb = 0.0,
    cykerlb = 0.0,
    cxkerub = 1.0,
    cykerub = 1.0
  )

  fast.u <- getFromNamespace(".np_inid_boot_from_ksum_unconditional", "np")
  fast.c <- getFromNamespace(".np_inid_boot_from_ksum_conditional", "np")

  u.dens.fast <- fast.u(
    xdat = tx,
    exdat = ex,
    bws = u.dens.bw,
    B = B,
    operator = "normal",
    counts = counts
  )
  u.dist.fast <- fast.u(
    xdat = tx,
    exdat = ex,
    bws = u.dist.bw,
    B = B,
    operator = "integral",
    counts = counts
  )
  c.dens.fast <- fast.c(
    xdat = tx,
    ydat = ty,
    exdat = ex,
    eydat = ey,
    bws = c.dens.bw,
    B = B,
    cdf = FALSE,
    counts = counts
  )
  c.dist.fast <- fast.c(
    xdat = tx,
    ydat = ty,
    exdat = ex,
    eydat = ey,
    bws = c.dist.bw,
    B = B,
    cdf = TRUE,
    counts = counts
  )

  u.dens.explicit <- matrix(NA_real_, nrow = B, ncol = nrow(ex))
  u.dist.explicit <- matrix(NA_real_, nrow = B, ncol = nrow(ex))
  c.dens.explicit <- matrix(NA_real_, nrow = B, ncol = nrow(ex))
  c.dist.explicit <- matrix(NA_real_, nrow = B, ncol = nrow(ex))
  for (b in seq_len(B)) {
    idx <- rep.int(seq_len(n), counts[, b])
    u.dens.explicit[b, ] <- npudens(tdat = tx[idx, , drop = FALSE], edat = ex, bws = u.dens.bw)$dens
    u.dist.explicit[b, ] <- npudist(tdat = tx[idx, , drop = FALSE], edat = ex, bws = u.dist.bw)$dist
    c.dens.explicit[b, ] <- npcdens(
      txdat = tx[idx, , drop = FALSE],
      tydat = ty[idx, , drop = FALSE],
      exdat = ex,
      eydat = ey,
      bws = c.dens.bw
    )$condens
    c.dist.explicit[b, ] <- npcdist(
      txdat = tx[idx, , drop = FALSE],
      tydat = ty[idx, , drop = FALSE],
      exdat = ex,
      eydat = ey,
      bws = c.dist.bw
    )$condist
  }

  expect_equal(u.dens.fast$t, u.dens.explicit, tolerance = 1e-10)
  expect_equal(u.dist.fast$t, u.dist.explicit, tolerance = 1e-10)
  expect_equal(c.dens.fast$t, c.dens.explicit, tolerance = 1e-10)
  expect_equal(c.dist.fast$t, c.dist.explicit, tolerance = 1e-10)
})

test_that("conditional helper constructors forward kernel and bwscaling options", {
  skip_if_not_installed("np")

  set.seed(3262)
  n <- 55
  tx <- data.frame(x = runif(n))
  ty <- data.frame(y = runif(n))
  bw.base <- npcdensbw(
    xdat = tx,
    ydat = ty,
    bws = c(0.22, 0.22),
    bandwidth.compute = FALSE,
    bwtype = "fixed",
    bwscaling = FALSE,
    cxkertype = "epanechnikov",
    cykertype = "epanechnikov",
    cxkerbound = "fixed",
    cykerbound = "fixed",
    cxkerlb = 0.0,
    cykerlb = 0.0,
    cxkerub = 1.0,
    cykerub = 1.0
  )
  bw <- bw.base
  bw$bwscaling <- TRUE

  np.ns <- asNamespace("np")
  cap <- new.env(parent = emptyenv())
  cap$calls <- list()

  trace(
    what = "kbandwidth.numeric",
    where = np.ns,
    tracer = bquote({
      assign(
        "calls",
        c(
          get("calls", envir = .(cap)),
          list(list(
            bwscaling = bwscaling,
            ckertype = ckertype,
            ckerorder = ckerorder,
            ckerbound = ckerbound,
            ckerlb = ckerlb,
            ckerub = ckerub
          ))
        ),
        envir = .(cap)
      )
    }),
    print = FALSE
  )
  on.exit(untrace("kbandwidth.numeric", where = np.ns), add = TRUE)

  make.kx <- getFromNamespace(".np_con_make_kbandwidth_x", "np")
  make.kxy <- getFromNamespace(".np_con_make_kbandwidth_xy", "np")
  kx <- make.kx(bws = bw, xdat = tx)
  kxy <- make.kxy(bws = bw, xdat = tx, ydat = ty)
  expect_false(is.null(kx))
  expect_false(is.null(kxy))
  expect_true(length(cap$calls) >= 2L)

  for (call in cap$calls) {
    expect_identical(isTRUE(call$bwscaling), FALSE)
    expect_identical(as.character(call$ckertype), as.character(bw$cxkertype))
    expect_identical(as.character(call$ckerorder), as.character(bw$cxkerorder))
    expect_identical(as.character(call$ckerbound), as.character(bw$cxkerbound))
  }
})

test_that("manual bws bwscaling toggle is invariant for unconditional ksum helper output", {
  skip_if_not_installed("np")

  set.seed(3263)
  n <- 70
  tx <- data.frame(x = runif(n))
  ex <- data.frame(x = seq(0.04, 0.96, length.out = 11))
  B <- 7L
  counts <- rmultinom(n = B, size = n, prob = rep.int(1 / n, n))

  bw.base <- npudensbw(
    dat = tx,
    bws = 0.21,
    bandwidth.compute = FALSE,
    bwtype = "fixed",
    bwscaling = FALSE
  )
  bw.alt <- bw.base
  bw.alt$bwscaling <- TRUE

  fast.u <- getFromNamespace(".np_inid_boot_from_ksum_unconditional", "np")
  out.base <- fast.u(
    xdat = tx,
    exdat = ex,
    bws = bw.base,
    B = B,
    operator = "normal",
    counts = counts
  )
  out.alt <- fast.u(
    xdat = tx,
    exdat = ex,
    bws = bw.alt,
    B = B,
    operator = "normal",
    counts = counts
  )

  expect_equal(out.base$t0, out.alt$t0, tolerance = 1e-12)
  expect_equal(out.base$t, out.alt$t, tolerance = 1e-12)
})

test_that("lp degree is used by regression inid helper construction", {
  skip_if_not_installed("np")

  set.seed(3264)
  n <- 75
  x1 <- runif(n)
  x2 <- runif(n)
  tx <- data.frame(x1 = x1, x2 = x2)
  y <- sin(3 * x1) + 0.25 * x2^2 + rnorm(n, sd = 0.04)
  ex <- tx[seq_len(20), , drop = FALSE]
  B <- 6L
  counts <- rmultinom(n = B, size = n, prob = rep.int(1 / n, n))

  bw2 <- npregbw(
    xdat = tx,
    ydat = y,
    regtype = "lp",
    basis = "glp",
    degree = c(2L, 2L),
    bws = c(0.28, 0.28),
    bandwidth.compute = FALSE
  )
  bw3 <- npregbw(
    xdat = tx,
    ydat = y,
    regtype = "lp",
    basis = "glp",
    degree = c(3L, 3L),
    bws = c(0.28, 0.28),
    bandwidth.compute = FALSE
  )

  fast.fun <- getFromNamespace(".np_inid_boot_from_regression", "np")
  out2 <- fast.fun(xdat = tx, exdat = ex, bws = bw2, ydat = y, B = B, counts = counts)
  out3 <- fast.fun(xdat = tx, exdat = ex, bws = bw3, ydat = y, B = B, counts = counts)

  fit2 <- npreg(txdat = tx, tydat = y, exdat = ex, bws = bw2, gradients = FALSE, warn.glp.gradient = FALSE)$mean
  fit3 <- npreg(txdat = tx, tydat = y, exdat = ex, bws = bw3, gradients = FALSE, warn.glp.gradient = FALSE)$mean

  expect_equal(out2$t0, fit2, tolerance = 1e-6)
  expect_equal(out3$t0, fit3, tolerance = 1e-6)
  expect_gt(max(abs(out2$t0 - out3$t0)), 1e-6)
})

test_that("npreg plot bootstrap inid supports lp basis variants", {
  skip_if_not_installed("np")

  set.seed(327)
  n <- 70
  x <- runif(n)
  y <- sin(2 * pi * x) + rnorm(n, sd = 0.1)
  tx <- data.frame(x = x)

  for (basis in c("additive", "tensor", "glp")) {
    bw <- npregbw(
      xdat = tx,
      ydat = y,
      regtype = "lp",
      degree = 2,
      basis = basis,
      bws = 0.2,
      bandwidth.compute = FALSE
    )

    out <- suppressWarnings(
      plot(
        bw,
        xdat = tx,
        ydat = y,
        plot.behavior = "data",
        perspective = FALSE,
        plot.errors.method = "bootstrap",
        plot.errors.boot.method = "inid",
        plot.errors.boot.num = 7
      )
    )
    expect_type(out, "list")
    expect_true(length(out) > 0)
  }
})

test_that("density/distribution plot bootstrap rejects wild selector", {
  skip_if_not_installed("np")

  set.seed(328)
  n <- 50
  x <- rnorm(n)
  y <- rnorm(n)

  ubw <- npudensbw(dat = data.frame(x = x), bws = 0.5, bandwidth.compute = FALSE)
  expect_error(
    suppressWarnings(
      plot(
        ubw,
        plot.behavior = "data",
        perspective = FALSE,
        plot.errors.method = "bootstrap",
        plot.errors.boot.method = "wild",
        plot.errors.boot.num = 9
      )
    ),
    "not supported for unconditional density/distribution estimators"
  )

  cbw <- npcdensbw(
    xdat = data.frame(x = x),
    ydat = data.frame(y = y),
    bws = c(0.6, 0.6),
    bandwidth.compute = FALSE
  )
  expect_error(
    suppressWarnings(
      plot(
        cbw,
        plot.behavior = "data",
        perspective = FALSE,
        plot.errors.method = "bootstrap",
        plot.errors.boot.method = "wild",
        plot.errors.boot.num = 9
      )
    ),
    "not supported for conditional density/distribution estimators"
  )
})

test_that("nearest-neighbor plot helpers run for regression and density/distribution families", {
  skip_if_not_installed("np")

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

test_that("nearest-neighbor frozen bootstrap plots run for unsupervised families and reject unsupported regression families", {
  skip_if_not_installed("np")

  set.seed(32908)
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
      plot.errors.boot.nonfixed = "frozen",
      plot.errors.boot.num = 5,
      ...
    ))
  }

  for (bt in c("generalized_nn", "adaptive_nn")) {
    rbw <- npregbw(xdat = x, ydat = y, regtype = "ll", nmulti = 1, bwtype = bt)
    expect_error(
      run_plot(rbw, xdat = x, ydat = y, plot.errors.boot.method = "inid"),
      "currently supported only"
    )

    for (boot.method in c("inid", "fixed", "geom")) {
      ubw <- npudensbw(dat = x, nmulti = 1, bwtype = bt)
      expect_type(run_plot(ubw, plot.errors.boot.method = boot.method), "list")

      dbw <- npudistbw(dat = x, nmulti = 1, bwtype = bt)
      expect_type(run_plot(dbw, plot.errors.boot.method = boot.method), "list")

      cbw <- npcdensbw(xdat = x, ydat = yframe, nmulti = 1, bwtype = bt)
      expect_type(
        run_plot(cbw, xdat = x, ydat = yframe, view = "fixed", plot.errors.boot.method = boot.method),
        "list"
      )

      cdbw <- npcdistbw(xdat = x, ydat = yframe, nmulti = 1, bwtype = bt)
      expect_type(
        run_plot(cdbw, xdat = x, ydat = yframe, view = "fixed", plot.errors.boot.method = boot.method),
        "list"
      )
    }
  }
})

test_that("nearest-neighbor regression helper matches explicit resample refits", {
  skip_if_not_installed("np")

  set.seed(32905)
  n <- 45
  x1 <- runif(n)
  x2 <- runif(n)
  y <- sin(2 * pi * x1) + 0.6 * x2 + rnorm(n, sd = 0.08)
  tx <- data.frame(x1 = x1, x2 = x2)
  ex <- tx[seq_len(15), , drop = FALSE]
  B <- 9L
  counts <- rmultinom(n = B, size = n, prob = rep.int(1 / n, n))
  fast.fun <- getFromNamespace(".np_inid_boot_from_reghat_exact", "np")

  cfgs <- list(
    list(regtype = "lc", basis = NULL, degree = NULL),
    list(regtype = "ll", basis = NULL, degree = NULL),
    list(regtype = "lp", basis = "glp", degree = c(2L, 2L))
  )

  for (bt in c("generalized_nn", "adaptive_nn")) {
    bw.val <- if (identical(bt, "adaptive_nn")) c(5, 5) else c(1, 1)

    for (cfg in cfgs) {
      bw.args <- list(
        xdat = tx,
        ydat = y,
        regtype = cfg$regtype,
        bws = bw.val,
        bwtype = bt,
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

      expect_equal(
        fast.out$t,
        explicit.t,
        tolerance = 1e-10,
        info = paste(bt, cfg$regtype)
      )
      expect_equal(
        fast.out$t0,
        npreg(txdat = tx, tydat = y, exdat = ex, bws = bw, gradients = FALSE, warn.glp.gradient = FALSE)$mean,
        tolerance = 1e-10,
        info = paste(bt, cfg$regtype)
      )
    }

    bw.grad <- npregbw(
      xdat = tx,
      ydat = y,
      regtype = "ll",
      bws = bw.val,
      bwtype = bt,
      bandwidth.compute = FALSE
    )
    fast.grad <- fast.fun(
      xdat = tx,
      exdat = ex,
      bws = bw.grad,
      ydat = y,
      B = B,
      counts = counts,
      gradients = TRUE,
      gradient.order = 1L,
      slice.index = 1L
    )

    explicit.grad <- matrix(NA_real_, nrow = B, ncol = nrow(ex))
    for (b in seq_len(B)) {
      idx <- rep.int(seq_len(n), counts[, b])
      explicit.grad[b, ] <- npreg(
        txdat = tx[idx, , drop = FALSE],
        tydat = y[idx],
        exdat = ex,
        bws = bw.grad,
        gradients = TRUE,
        gradient.order = 1L,
        warn.glp.gradient = FALSE
      )$grad[, 1L]
    }

    expect_equal(
      fast.grad$t,
      explicit.grad,
      tolerance = 1e-10,
      info = paste(bt, "ll-grad")
    )
    expect_equal(
      fast.grad$t0,
      npreg(txdat = tx, tydat = y, exdat = ex, bws = bw.grad, gradients = TRUE, gradient.order = 1L, warn.glp.gradient = FALSE)$grad[, 1L],
      tolerance = 1e-10,
      info = paste(bt, "ll-grad")
    )
  }
})

test_that("nearest-neighbor ksum helpers match explicit resample refits", {
  skip_if_not_installed("np")

  set.seed(32906)
  n <- 40
  tx <- data.frame(x = runif(n))
  ty <- data.frame(y = runif(n))
  ex <- data.frame(x = seq(0.1, 0.9, length.out = 11))
  ey <- data.frame(y = seq(0.1, 0.9, length.out = 11))
  B <- 9L
  counts <- rmultinom(n = B, size = n, prob = rep.int(1 / n, n))

  fast.u <- getFromNamespace(".np_inid_boot_from_ksum_unconditional", "np")
  fast.c <- getFromNamespace(".np_inid_boot_from_ksum_conditional", "np")

  for (bt in c("generalized_nn", "adaptive_nn")) {
    bw.val <- if (identical(bt, "adaptive_nn")) 5 else 1
    cbw.val <- rep.int(bw.val, 2L)

    u.dens.bw <- npudensbw(dat = tx, bws = bw.val, bwtype = bt, bandwidth.compute = FALSE)
    u.dist.bw <- npudistbw(dat = tx, bws = bw.val, bwtype = bt, bandwidth.compute = FALSE)
    c.dens.bw <- npcdensbw(xdat = tx, ydat = ty, bws = cbw.val, bwtype = bt, bandwidth.compute = FALSE)
    c.dist.bw <- npcdistbw(xdat = tx, ydat = ty, bws = cbw.val, bwtype = bt, bandwidth.compute = FALSE)

    u.dens.fast <- fast.u(
      xdat = tx,
      exdat = ex,
      bws = u.dens.bw,
      B = B,
      operator = "normal",
      counts = counts
    )
    u.dist.fast <- fast.u(
      xdat = tx,
      exdat = ex,
      bws = u.dist.bw,
      B = B,
      operator = "integral",
      counts = counts
    )
    c.dens.fast <- fast.c(
      xdat = tx,
      ydat = ty,
      exdat = ex,
      eydat = ey,
      bws = c.dens.bw,
      B = B,
      cdf = FALSE,
      counts = counts
    )
    c.dist.fast <- fast.c(
      xdat = tx,
      ydat = ty,
      exdat = ex,
      eydat = ey,
      bws = c.dist.bw,
      B = B,
      cdf = TRUE,
      counts = counts
    )

    u.dens.explicit <- matrix(NA_real_, nrow = B, ncol = nrow(ex))
    u.dist.explicit <- matrix(NA_real_, nrow = B, ncol = nrow(ex))
    c.dens.explicit <- matrix(NA_real_, nrow = B, ncol = nrow(ex))
    c.dist.explicit <- matrix(NA_real_, nrow = B, ncol = nrow(ex))
    for (b in seq_len(B)) {
      idx <- rep.int(seq_len(n), counts[, b])
      u.dens.explicit[b, ] <- npudens(tdat = tx[idx, , drop = FALSE], edat = ex, bws = u.dens.bw)$dens
      u.dist.explicit[b, ] <- npudist(tdat = tx[idx, , drop = FALSE], edat = ex, bws = u.dist.bw)$dist
      c.dens.explicit[b, ] <- npcdens(
        txdat = tx[idx, , drop = FALSE],
        tydat = ty[idx, , drop = FALSE],
        exdat = ex,
        eydat = ey,
        bws = c.dens.bw
      )$condens
      c.dist.explicit[b, ] <- npcdist(
        txdat = tx[idx, , drop = FALSE],
        tydat = ty[idx, , drop = FALSE],
        exdat = ex,
        eydat = ey,
        bws = c.dist.bw
      )$condist
    }

    expect_equal(u.dens.fast$t, u.dens.explicit, tolerance = 1e-10, info = paste(bt, "udens"))
    expect_equal(u.dist.fast$t, u.dist.explicit, tolerance = 1e-10, info = paste(bt, "udist"))
    expect_equal(c.dens.fast$t, c.dens.explicit, tolerance = 1e-10, info = paste(bt, "npcdens"))
    expect_equal(c.dist.fast$t, c.dist.explicit, tolerance = 1e-10, info = paste(bt, "npcdist"))
    expect_equal(u.dens.fast$t0, npudens(tdat = tx, edat = ex, bws = u.dens.bw)$dens, tolerance = 1e-12, info = paste(bt, "udens-t0"))
    expect_equal(u.dist.fast$t0, npudist(tdat = tx, edat = ex, bws = u.dist.bw)$dist, tolerance = 1e-12, info = paste(bt, "udist-t0"))
    expect_equal(c.dens.fast$t0, npcdens(txdat = tx, tydat = ty, exdat = ex, eydat = ey, bws = c.dens.bw)$condens, tolerance = 1e-12, info = paste(bt, "npcdens-t0"))
    expect_equal(c.dist.fast$t0, npcdist(txdat = tx, tydat = ty, exdat = ex, eydat = ey, bws = c.dist.bw)$condist, tolerance = 1e-12, info = paste(bt, "npcdist-t0"))
  }
})

test_that("nearest-neighbor frozen hat helpers match frozen operator application", {
  skip_if_not_installed("np")

  set.seed(32907)
  n <- 36
  tx <- data.frame(x = runif(n))
  ty <- data.frame(y = runif(n))
  ex <- data.frame(x = seq(0.1, 0.9, length.out = 9))
  ey <- data.frame(y = seq(0.15, 0.85, length.out = 9))
  B <- 7L
  counts <- rmultinom(n = B, size = n, prob = rep.int(1 / n, n))

  frozen.u <- getFromNamespace(".np_inid_boot_from_hat_unconditional_frozen", "np")
  frozen.c <- getFromNamespace(".np_inid_boot_from_hat_conditional_frozen", "np")
  udens.hat <- getFromNamespace("npudenshat", "np")
  udist.hat <- getFromNamespace("npudisthat", "np")
  cdens.hat <- getFromNamespace("npcdenshat", "np")
  cdist.hat <- getFromNamespace("npcdisthat", "np")

  for (bt in c("generalized_nn", "adaptive_nn")) {
    bw.val <- if (identical(bt, "adaptive_nn")) 5 else 1
    cbw.val <- rep.int(bw.val, 2L)

    u.dens.bw <- npudensbw(dat = tx, bws = bw.val, bwtype = bt, bandwidth.compute = FALSE)
    u.dist.bw <- npudistbw(dat = tx, bws = bw.val, bwtype = bt, bandwidth.compute = FALSE)
    c.dens.bw <- npcdensbw(xdat = tx, ydat = ty, bws = cbw.val, bwtype = bt, bandwidth.compute = FALSE)
    c.dist.bw <- npcdistbw(xdat = tx, ydat = ty, bws = cbw.val, bwtype = bt, bandwidth.compute = FALSE)

    H.udens <- unclass(udens.hat(bws = u.dens.bw, tdat = tx, edat = ex, output = "matrix"))
    H.udist <- unclass(udist.hat(bws = u.dist.bw, tdat = tx, edat = ex, output = "matrix"))
    H.cdens <- unclass(cdens.hat(bws = c.dens.bw, txdat = tx, tydat = ty, exdat = ex, eydat = ey, output = "matrix"))
    H.cdist <- unclass(cdist.hat(bws = c.dist.bw, txdat = tx, tydat = ty, exdat = ex, eydat = ey, output = "matrix"))

    u.dens.frozen <- frozen.u(
      xdat = tx,
      exdat = ex,
      bws = u.dens.bw,
      B = B,
      operator = "normal",
      counts = counts
    )
    u.dist.frozen <- frozen.u(
      xdat = tx,
      exdat = ex,
      bws = u.dist.bw,
      B = B,
      operator = "integral",
      counts = counts
    )
    c.dens.frozen <- frozen.c(
      xdat = tx,
      ydat = ty,
      exdat = ex,
      eydat = ey,
      bws = c.dens.bw,
      B = B,
      cdf = FALSE,
      counts = counts
    )
    c.dist.frozen <- frozen.c(
      xdat = tx,
      ydat = ty,
      exdat = ex,
      eydat = ey,
      bws = c.dist.bw,
      B = B,
      cdf = TRUE,
      counts = counts
    )

    expect_equal(u.dens.frozen$t, t(H.udens %*% counts), tolerance = 1e-10, info = paste(bt, "udens"))
    expect_equal(u.dist.frozen$t, t(H.udist %*% counts), tolerance = 1e-10, info = paste(bt, "udist"))
    expect_equal(c.dens.frozen$t, t(H.cdens %*% counts), tolerance = 1e-10, info = paste(bt, "npcdens"))
    expect_equal(c.dist.frozen$t, t(H.cdist %*% counts), tolerance = 1e-10, info = paste(bt, "npcdist"))

    expect_equal(u.dens.frozen$t0, rowSums(H.udens), tolerance = 1e-12, info = paste(bt, "udens-t0"))
    expect_equal(u.dist.frozen$t0, rowSums(H.udist), tolerance = 1e-12, info = paste(bt, "udist-t0"))
    expect_equal(c.dens.frozen$t0, rowSums(H.cdens), tolerance = 1e-12, info = paste(bt, "npcdens-t0"))
    expect_equal(c.dist.frozen$t0, rowSums(H.cdist), tolerance = 1e-12, info = paste(bt, "npcdist-t0"))
  }
})

test_that("npindex plot bootstrap inid supports nearest-neighbor bwtypes", {
  skip_if_not_installed("np")

  set.seed(3291)
  n <- 40
  x1 <- runif(n)
  x2 <- runif(n)
  y <- sin(x1 + x2) + rnorm(n, sd = 0.08)
  tx <- data.frame(x1 = x1, x2 = x2)

  run_plot <- function(bw) {
    suppressWarnings(plot(
      bw,
      xdat = tx,
      ydat = y,
      plot.behavior = "data",
      perspective = FALSE,
      plot.errors.method = "bootstrap",
      plot.errors.boot.method = "inid",
      plot.errors.boot.num = 5
    ))
  }

  for (bt in c("generalized_nn", "adaptive_nn")) {
    bw.lc <- npindexbw(
      xdat = tx,
      ydat = y,
      bws = c(1, 1, 5L),
      bandwidth.compute = FALSE,
      bwtype = bt,
      regtype = "lc"
    )
    out.lc <- run_plot(bw.lc)
    expect_type(out.lc, "list")
    expect_true(length(out.lc) > 0)

    bw.ll <- npindexbw(
      xdat = tx,
      ydat = y,
      bws = c(1, 1, 5L),
      bandwidth.compute = FALSE,
      bwtype = bt,
      regtype = "ll"
    )
    out.ll <- run_plot(bw.ll)
    expect_type(out.ll, "list")
    expect_true(length(out.ll) > 0)
  }
})

test_that("npindex nearest-neighbor inid helper matches explicit resample refits", {
  skip_if_not_installed("np")

  set.seed(3292)
  n <- 45
  x1 <- runif(n)
  x2 <- runif(n)
  y <- sin(x1 + x2) + rnorm(n, sd = 0.08)
  tx <- data.frame(x1 = x1, x2 = x2)
  B <- 9L
  counts <- rmultinom(n = B, size = n, prob = rep.int(1 / n, n))
  fast.fun <- getFromNamespace(".np_inid_boot_from_index", "np")

  cfgs <- list(
    list(regtype = "lc", basis = NULL, degree = NULL, h = 5L),
    list(regtype = "ll", basis = NULL, degree = NULL, h = 5L),
    list(regtype = "lp", basis = "tensor", degree = 2L, h = 5L)
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

      fast.out <- fast.fun(
        xdat = tx,
        ydat = y,
        bws = bw,
        B = B,
        counts = counts,
        gradients = FALSE
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

      expect_equal(
        fast.out$t,
        explicit.t,
        tolerance = 1e-8,
        info = paste(bt, cfg$regtype)
      )
      expect_equal(
        fast.out$t0,
        npindex(txdat = tx, tydat = y, exdat = tx, bws = bw, gradients = FALSE)$mean,
        tolerance = 1e-8,
        info = paste(bt, cfg$regtype)
      )
    }
  }
})
