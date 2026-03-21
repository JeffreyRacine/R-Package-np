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
