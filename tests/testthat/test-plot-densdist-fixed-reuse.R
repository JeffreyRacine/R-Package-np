test_that("fixed unconditional density/distribution helper matches explicit refits", {
  skip_if_not_installed("np")

  set.seed(603101)
  n <- 70
  tx <- data.frame(x = rnorm(n))
  ex <- data.frame(x = seq(min(tx$x), max(tx$x), length.out = 19))
  B <- 7L
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

test_that("fixed conditional density/distribution helper matches explicit refits", {
  skip_if_not_installed("np")

  set.seed(603102)
  n <- 80
  tx <- data.frame(x = rnorm(n))
  ty <- data.frame(y = rnorm(n))
  ex <- data.frame(x = seq(min(tx$x), max(tx$x), length.out = 17))
  ey <- data.frame(y = seq(min(ty$y), max(ty$y), length.out = 17))
  B <- 7L
  counts <- rmultinom(n = B, size = n, prob = rep.int(1 / n, n))

  fast.fun <- getFromNamespace(".np_inid_boot_from_ksum_conditional", "np")
  oracle.fun <- getFromNamespace(".np_inid_boot_from_conditional_localpoly_fixed_rowwise", "np")

  cfgs <- list(
    list(name = "lc", regtype = "lc", degree = NULL),
    list(name = "ll", regtype = "ll", degree = NULL),
    list(name = "lp", regtype = "lp", degree = 2L)
  )

  for (cfg in cfgs) {
    dens.args <- list(
      xdat = tx,
      ydat = ty,
      bws = c(0.45, 0.45),
      bandwidth.compute = FALSE,
      bwtype = "fixed",
      regtype = cfg$regtype
    )
    dist.args <- dens.args
    if (!is.null(cfg$degree)) {
      dens.args$degree <- cfg$degree
      dens.args$basis <- "glp"
      dens.args$bernstein.basis <- FALSE
      dist.args$degree <- cfg$degree
      dist.args$basis <- "glp"
      dist.args$bernstein.basis <- FALSE
    }

    dens.bw <- do.call(npcdensbw, dens.args)
    dist.bw <- do.call(npcdistbw, dist.args)
    oracle.t0.tolerance <- if (identical(cfg$regtype, "lp")) 1e-10 else 1e-12

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

    expect_equal(
      dens.fast$t,
      dens.explicit,
      tolerance = 1e-10,
      info = paste("density explicit", cfg$name)
    )
    expect_equal(
      dist.fast$t,
      dist.explicit,
      tolerance = 1e-10,
      info = paste("distribution explicit", cfg$name)
    )
    expect_equal(
      dens.fast$t0,
      npcdens(txdat = tx, tydat = ty, exdat = ex, eydat = ey, bws = dens.bw)$condens,
      tolerance = 1e-12,
      info = paste("density t0", cfg$name)
    )
    expect_equal(
      dist.fast$t0,
      npcdist(txdat = tx, tydat = ty, exdat = ex, eydat = ey, bws = dist.bw)$condist,
      tolerance = 1e-12,
      info = paste("distribution t0", cfg$name)
    )

    if (!identical(cfg$regtype, "lc")) {
      dens.oracle <- oracle.fun(
        xdat = tx,
        ydat = ty,
        exdat = ex,
        eydat = ey,
        bws = dens.bw,
        B = B,
        cdf = FALSE,
        counts = counts
      )
      dist.oracle <- oracle.fun(
        xdat = tx,
        ydat = ty,
        exdat = ex,
        eydat = ey,
        bws = dist.bw,
        B = B,
        cdf = TRUE,
        counts = counts
      )

      expect_equal(
        dens.fast$t,
        dens.oracle$t,
        tolerance = 1e-12,
        info = paste("density oracle", cfg$name)
      )
      expect_equal(
        dist.fast$t,
        dist.oracle$t,
        tolerance = 1e-12,
        info = paste("distribution oracle", cfg$name)
      )
    expect_equal(
      dens.fast$t0,
      dens.oracle$t0,
      tolerance = oracle.t0.tolerance,
      info = paste("density oracle t0", cfg$name)
    )
    expect_equal(
      dist.fast$t0,
      dist.oracle$t0,
      tolerance = oracle.t0.tolerance,
      info = paste("distribution oracle t0", cfg$name)
    )

      drawer <- function(start, stopi) counts[, start:stopi, drop = FALSE]
      dens.fast.drawer <- fast.fun(
        xdat = tx,
        ydat = ty,
        exdat = ex,
        eydat = ey,
        bws = dens.bw,
        B = B,
        cdf = FALSE,
        counts.drawer = drawer
      )
      dist.fast.drawer <- fast.fun(
        xdat = tx,
        ydat = ty,
        exdat = ex,
        eydat = ey,
        bws = dist.bw,
        B = B,
        cdf = TRUE,
        counts.drawer = drawer
      )
      dens.oracle.drawer <- oracle.fun(
        xdat = tx,
        ydat = ty,
        exdat = ex,
        eydat = ey,
        bws = dens.bw,
        B = B,
        cdf = FALSE,
        counts.drawer = drawer
      )
      dist.oracle.drawer <- oracle.fun(
        xdat = tx,
        ydat = ty,
        exdat = ex,
        eydat = ey,
        bws = dist.bw,
        B = B,
        cdf = TRUE,
        counts.drawer = drawer
      )

      expect_equal(
        dens.fast.drawer$t,
        dens.oracle.drawer$t,
        tolerance = 1e-12,
        info = paste("density drawer oracle", cfg$name)
      )
      expect_equal(
        dist.fast.drawer$t,
        dist.oracle.drawer$t,
        tolerance = 1e-12,
        info = paste("distribution drawer oracle", cfg$name)
      )
    expect_equal(
      dens.fast.drawer$t0,
      dens.oracle.drawer$t0,
      tolerance = oracle.t0.tolerance,
      info = paste("density drawer oracle t0", cfg$name)
    )
    expect_equal(
      dist.fast.drawer$t0,
      dist.oracle.drawer$t0,
      tolerance = oracle.t0.tolerance,
      info = paste("distribution drawer oracle t0", cfg$name)
    )
    }
  }
})

test_that("fixed conditional ll/lp grouped inid helper matches the rowwise oracle on repeated evaluation x rows", {
  skip_if_not_installed("np")

  set.seed(6031021)
  n <- 72
  tx <- data.frame(x = rnorm(n))
  ty <- data.frame(y = rnorm(n))
  x.grid <- seq(min(tx$x), max(tx$x), length.out = 6L)
  y.grid <- seq(min(ty$y), max(ty$y), length.out = 5L)
  grid <- expand.grid(y = y.grid, x = x.grid)
  ex <- data.frame(x = grid$x)
  ey <- data.frame(y = grid$y)
  B <- 6L
  counts <- rmultinom(n = B, size = n, prob = rep.int(1 / n, n))

  fast.fun <- getFromNamespace(".np_inid_boot_from_ksum_conditional", "np")
  oracle.fun <- getFromNamespace(".np_inid_boot_from_conditional_localpoly_fixed_rowwise", "np")

  cfgs <- list(
    list(name = "ll", regtype = "ll", degree = NULL),
    list(name = "lp", regtype = "lp", degree = 2L)
  )

  for (cfg in cfgs) {
    dens.args <- list(
      xdat = tx,
      ydat = ty,
      bws = c(0.40, 0.40),
      bandwidth.compute = FALSE,
      bwtype = "fixed",
      regtype = cfg$regtype
    )
    dist.args <- dens.args
    if (!is.null(cfg$degree)) {
      dens.args$degree <- cfg$degree
      dens.args$basis <- "glp"
      dens.args$bernstein.basis <- FALSE
      dist.args$degree <- cfg$degree
      dist.args$basis <- "glp"
      dist.args$bernstein.basis <- FALSE
    }

    dens.bw <- do.call(npcdensbw, dens.args)
    dist.bw <- do.call(npcdistbw, dist.args)

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
    dens.oracle <- oracle.fun(
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
    dist.oracle <- oracle.fun(
      xdat = tx,
      ydat = ty,
      exdat = ex,
      eydat = ey,
      bws = dist.bw,
      B = B,
      cdf = TRUE,
      counts = counts
    )

    expect_equal(dens.fast$t, dens.oracle$t, tolerance = 1e-12, info = paste("density repeated-x oracle", cfg$name))
    expect_equal(dens.fast$t0, dens.oracle$t0, tolerance = 1e-10, info = paste("density repeated-x oracle t0", cfg$name))
    expect_equal(dist.fast$t, dist.oracle$t, tolerance = 1e-12, info = paste("distribution repeated-x oracle", cfg$name))
    expect_equal(dist.fast$t0, dist.oracle$t0, tolerance = 1e-10, info = paste("distribution repeated-x oracle t0", cfg$name))
  }
})

test_that("fixed conditional ll/lp grouped helper is drawer-equal and chunk-invariant on repeated evaluation x rows", {
  skip_if_not_installed("np")

  old.chunk <- getOption("np.plot.inid.chunk.size")
  options(np.plot.inid.chunk.size = 2L)
  on.exit(options(np.plot.inid.chunk.size = old.chunk), add = TRUE)

  set.seed(6031022)
  n <- 68
  tx <- data.frame(x = rnorm(n))
  ty <- data.frame(y = rnorm(n))
  x.grid <- seq(min(tx$x), max(tx$x), length.out = 7L)
  y.grid <- seq(min(ty$y), max(ty$y), length.out = 6L)
  grid <- expand.grid(y = y.grid, x = x.grid)
  ex <- data.frame(x = grid$x)
  ey <- data.frame(y = grid$y)
  B <- 7L
  counts <- rmultinom(n = B, size = n, prob = rep.int(1 / n, n))

  fast.fun <- getFromNamespace(".np_inid_boot_from_ksum_conditional", "np")
  oracle.fun <- getFromNamespace(".np_inid_boot_from_conditional_localpoly_fixed_rowwise", "np")

  cfgs <- list(
    list(name = "ll", regtype = "ll", degree = NULL),
    list(name = "lp", regtype = "lp", degree = 2L)
  )

  for (cfg in cfgs) {
    dens.args <- list(
      xdat = tx,
      ydat = ty,
      bws = c(0.38, 0.38),
      bandwidth.compute = FALSE,
      bwtype = "fixed",
      regtype = cfg$regtype
    )
    dist.args <- dens.args
    if (!is.null(cfg$degree)) {
      dens.args$degree <- cfg$degree
      dens.args$basis <- "glp"
      dens.args$bernstein.basis <- FALSE
      dist.args$degree <- cfg$degree
      dist.args$basis <- "glp"
      dist.args$bernstein.basis <- FALSE
    }

    dens.bw <- do.call(npcdensbw, dens.args)
    dist.bw <- do.call(npcdistbw, dist.args)
    drawer <- function(start, stopi) counts[, start:stopi, drop = FALSE]

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
    dens.drawer <- fast.fun(
      xdat = tx,
      ydat = ty,
      exdat = ex,
      eydat = ey,
      bws = dens.bw,
      B = B,
      cdf = FALSE,
      counts.drawer = drawer
    )
    dens.oracle.drawer <- oracle.fun(
      xdat = tx,
      ydat = ty,
      exdat = ex,
      eydat = ey,
      bws = dens.bw,
      B = B,
      cdf = FALSE,
      counts.drawer = drawer
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
    dist.drawer <- fast.fun(
      xdat = tx,
      ydat = ty,
      exdat = ex,
      eydat = ey,
      bws = dist.bw,
      B = B,
      cdf = TRUE,
      counts.drawer = drawer
    )
    dist.oracle.drawer <- oracle.fun(
      xdat = tx,
      ydat = ty,
      exdat = ex,
      eydat = ey,
      bws = dist.bw,
      B = B,
      cdf = TRUE,
      counts.drawer = drawer
    )

    expect_equal(dens.drawer$t, dens.fast$t, tolerance = 1e-12, info = paste("density drawer chunk invariance", cfg$name))
    expect_equal(dens.drawer$t0, dens.fast$t0, tolerance = 1e-10, info = paste("density drawer chunk invariance t0", cfg$name))
    expect_equal(dens.drawer$t, dens.oracle.drawer$t, tolerance = 1e-8, info = paste("density drawer oracle repeated-x", cfg$name))
    expect_equal(dens.drawer$t0, dens.oracle.drawer$t0, tolerance = 1e-10, info = paste("density drawer oracle repeated-x t0", cfg$name))

    expect_equal(dist.drawer$t, dist.fast$t, tolerance = 1e-12, info = paste("distribution drawer chunk invariance", cfg$name))
    expect_equal(dist.drawer$t0, dist.fast$t0, tolerance = 1e-10, info = paste("distribution drawer chunk invariance t0", cfg$name))
    expect_equal(dist.drawer$t, dist.oracle.drawer$t, tolerance = 1e-8, info = paste("distribution drawer oracle repeated-x", cfg$name))
    expect_equal(dist.drawer$t0, dist.oracle.drawer$t0, tolerance = 1e-10, info = paste("distribution drawer oracle repeated-x t0", cfg$name))
  }
})

test_that("fixed density/distribution helpers preserve bounded kernel options", {
  skip_if_not_installed("np")

  set.seed(603103)
  n <- 60
  tx <- data.frame(x = runif(n))
  ty <- data.frame(y = runif(n))
  ex <- data.frame(x = seq(0.02, 0.98, length.out = 13))
  ey <- data.frame(y = seq(0.03, 0.97, length.out = 13))
  B <- 5L
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

test_that("fixed-bwtype unsupervised plot bootstrap covers inid fixed and geom", {
  skip_if_not_installed("np")

  set.seed(603104)
  n <- 48
  xdat <- data.frame(x = rnorm(n))
  ydat <- data.frame(y = rnorm(n))

  u.dens.bw <- npudensbw(dat = xdat, bws = 0.30, bandwidth.compute = FALSE, bwtype = "fixed")
  u.dist.bw <- npudistbw(dat = xdat, bws = 0.30, bandwidth.compute = FALSE, bwtype = "fixed")
  c.dens.bw <- npcdensbw(xdat = xdat, ydat = ydat, bws = c(0.35, 0.35), bandwidth.compute = FALSE, bwtype = "fixed")
  c.dist.bw <- npcdistbw(xdat = xdat, ydat = ydat, bws = c(0.35, 0.35), bandwidth.compute = FALSE, bwtype = "fixed")

  run_unsup_plot <- function(bw, ..., boot.method) {
    suppressWarnings(plot(
      bw,
      plot.behavior = "data",
      perspective = FALSE,
      plot.errors.method = "bootstrap",
      plot.errors.boot.method = boot.method,
      plot.errors.boot.blocklen = 3L,
      plot.errors.boot.num = 5L,
      plot.errors.type = "pointwise",
      neval = 11L,
      ...
    ))
  }

  for (boot.method in c("inid", "fixed", "geom")) {
    expect_type(run_unsup_plot(u.dens.bw, xdat = xdat, boot.method = boot.method), "list")
    expect_type(run_unsup_plot(u.dist.bw, xdat = xdat, boot.method = boot.method), "list")
    expect_type(run_unsup_plot(c.dens.bw, xdat = xdat, ydat = ydat, view = "fixed", boot.method = boot.method), "list")
    expect_type(run_unsup_plot(c.dist.bw, xdat = xdat, ydat = ydat, view = "fixed", boot.method = boot.method), "list")
  }
})
