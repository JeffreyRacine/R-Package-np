test_that("npscoefhat reproduces npscoef fitted values and supports matrix RHS", {
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
})

test_that("npscoef and npscoefhat support ll/lp basis variants", {
  set.seed(2469)
  n <- 95
  x <- runif(n)
  z <- runif(n)
  y <- (0.4 + x) * sin(2 * pi * z) + rnorm(n, sd = 0.04)
  tx <- data.frame(x = x)
  tz <- data.frame(z = z)
  ex <- data.frame(x = seq(min(x), max(x), length.out = 32))
  ez <- data.frame(z = seq(min(z), max(z), length.out = 32))

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
      bws = 0.18,
      bandwidth.compute = FALSE,
      regtype = cfg$regtype
    )
    if (!is.null(cfg$basis)) {
      bw.args$basis <- cfg$basis
      bw.args$degree <- 2L
    }
    bw <- do.call(npscoefbw, bw.args)
    expect_identical(bw$regtype, cfg$regtype, info = cfg$label)
    if (!is.null(cfg$basis))
      expect_identical(bw$basis, cfg$basis, info = cfg$label)

    fit.eval <- npscoef(
      bws = bw,
      txdat = tx,
      tydat = y,
      tzdat = tz,
      exdat = ex,
      ezdat = ez,
      errors = FALSE,
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
    expect_equal(as.vector(H.eval %*% y), as.vector(fit.eval$mean), tolerance = 1e-8, info = cfg$label)
  }
})

test_that("npscoefhat leave.one.out honors lc and ll paths", {
  set.seed(97533)
  n <- 85
  x <- runif(n)
  z <- runif(n)
  y <- (0.5 + x) * cos(2 * pi * z) + rnorm(n, sd = 0.04)
  tx <- data.frame(x = x)
  tz <- data.frame(z = z)
  ex <- data.frame(x = seq(min(x), max(x), length.out = 20))
  ez <- data.frame(z = seq(min(z), max(z), length.out = 20))

  bw.lc <- npscoefbw(
    xdat = x, zdat = z, ydat = y,
    bws = 0.16, regtype = "lc", bandwidth.compute = FALSE
  )
  H0.lc <- npscoefhat(
    bws = bw.lc, txdat = tx, tzdat = tz,
    output = "matrix", iterate = FALSE, leave.one.out = FALSE
  )
  H1.lc <- npscoefhat(
    bws = bw.lc, txdat = tx, tzdat = tz,
    output = "matrix", iterate = FALSE, leave.one.out = TRUE
  )
  expect_gt(max(abs(H0.lc - H1.lc)), 1e-8)

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

  expect_error(
    npscoefhat(
      bws = bw.ll,
      txdat = tx, tzdat = tz,
      exdat = ex, ezdat = ez,
      output = "matrix",
      iterate = FALSE,
      leave.one.out = TRUE
    ),
    "requires evaluation 'z' data to match training 'z' data"
  )
})

test_that("npplreghat reproduces npplreg fitted values and supports matrix RHS", {
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

test_that("npplreghat supports ll/lp with lp basis variants", {
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

test_that("semihat helper routes preserve LP and kernel-bound option contracts", {
  set.seed(314161)
  n <- 75
  x1 <- runif(n)
  x2 <- runif(n)
  z <- runif(n)
  y <- sin(2 * pi * z) + 0.8 * x1 - 0.3 * x2 + rnorm(n, sd = 0.04)

  tx <- data.frame(x1 = x1, x2 = x2)
  tz <- data.frame(z = z)
  ex <- tx[seq_len(25), , drop = FALSE]
  ez <- tz[seq_len(25), , drop = FALSE]

  sibw <- npindexbw(
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
    ckerbound = "none"
  )
  expect_identical(sibw$regtype, "lp")
  expect_identical(sibw$basis, "tensor")
  expect_true(isTRUE(sibw$bernstein.basis))
  expect_identical(sibw$ckerbound, "none")

  si.fit <- npindex(
    bws = sibw,
    txdat = tx,
    tydat = y,
    exdat = ex,
    gradients = FALSE
  )
  expect_true(all(is.finite(si.fit$mean)))
  si.H <- npindexhat(
    bws = sibw,
    txdat = tx,
    exdat = ex,
    s = 0L,
    output = "matrix"
  )
  si.apply <- npindexhat(
    bws = sibw,
    txdat = tx,
    exdat = ex,
    y = y,
    s = 0L,
    output = "apply"
  )
  expect_equal(as.vector(si.H %*% y), as.vector(si.apply), tolerance = 1e-10)

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
    "Invalid bounds for 'ckerbound'"
  )

  scbw <- npscoefbw(
    xdat = tx[, 1, drop = FALSE],
    zdat = tz,
    ydat = y,
    bws = 0.20,
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
  expect_identical(scbw$regtype, "lp")
  expect_identical(scbw$basis, "tensor")
  expect_true(isTRUE(scbw$bernstein.basis))
  expect_identical(scbw$ckerbound, "fixed")

  sc.fit <- npscoef(
    bws = scbw,
    txdat = tx[, 1, drop = FALSE],
    tydat = y,
    tzdat = tz,
    exdat = ex[, 1, drop = FALSE],
    ezdat = ez,
    iterate = FALSE
  )
  sc.H <- npscoefhat(
    bws = scbw,
    txdat = tx[, 1, drop = FALSE],
    tzdat = tz,
    exdat = ex[, 1, drop = FALSE],
    ezdat = ez,
    output = "matrix",
    iterate = FALSE
  )
  expect_equal(as.vector(sc.H %*% y), as.vector(sc.fit$mean), tolerance = 1e-8)

  plbw <- npplregbw(
    xdat = tx[, 1, drop = FALSE],
    zdat = tz,
    ydat = y,
    bws = matrix(c(0.20, 0.20), nrow = 2, ncol = 1),
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
  expect_identical(plbw$regtype, "lp")
  expect_identical(plbw$basis, "tensor")
  expect_true(isTRUE(plbw$bernstein.basis))
  expect_identical(plbw$ckerbound, "fixed")

  pl.fit <- npplreg(
    bws = plbw,
    txdat = tx[, 1, drop = FALSE],
    tydat = y,
    tzdat = tz,
    exdat = ex[, 1, drop = FALSE],
    ezdat = ez
  )
  pl.H <- npplreghat(
    bws = plbw,
    txdat = tx[, 1, drop = FALSE],
    tzdat = tz,
    exdat = ex[, 1, drop = FALSE],
    ezdat = ez,
    output = "matrix"
  )
  expect_true(all(is.finite(pl.fit$mean)))
  pl.apply <- npplreghat(
    bws = plbw,
    txdat = tx[, 1, drop = FALSE],
    tzdat = tz,
    exdat = ex[, 1, drop = FALSE],
    ezdat = ez,
    y = y,
    output = "apply"
  )
  expect_equal(as.vector(pl.H %*% y), as.vector(pl.apply), tolerance = 1e-10)
})

test_that("semihat validates class and scalar controls", {
  set.seed(27182)
  n <- 40
  x <- runif(n)
  y <- rnorm(n)
  z <- runif(n)

  rbw <- npregbw(y ~ x, bws = 0.2, bandwidth.compute = FALSE)
  expect_error(npindexhat(bws = rbw, txdat = data.frame(x = x)), "sibandwidth")
  expect_error(npplreghat(bws = rbw, txdat = data.frame(x = x), tzdat = data.frame(z = z)), "plbandwidth")
  expect_error(npscoefhat(bws = rbw, txdat = data.frame(x = x), tzdat = data.frame(z = z)), "scbandwidth")

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

  sc.cfgs <- list(
    list(regtype = "ll", basis = NULL, label = "ll"),
    list(regtype = "lp", basis = "glp", label = "lp-glp"),
    list(regtype = "lp", basis = "additive", label = "lp-additive"),
    list(regtype = "lp", basis = "tensor", label = "lp-tensor")
  )
  for (cfg in sc.cfgs) {
    sc.args <- list(
      xdat = x,
      zdat = z,
      ydat = y,
      bws = 0.2,
      bandwidth.compute = FALSE,
      regtype = cfg$regtype
    )
    if (!is.null(cfg$basis)) {
      sc.args$basis <- cfg$basis
      sc.args$degree <- 2L
    }
    scbw.cfg <- do.call(npscoefbw, sc.args)
    sc.out.cfg <- suppressWarnings(
      plot(
        scbw.cfg,
        xdat = data.frame(x = x),
        ydat = y,
        zdat = data.frame(z = z),
        perspective = FALSE,
        plot.behavior = "data",
        plot.errors.method = "bootstrap",
        plot.errors.boot.method = "wild",
        plot.errors.boot.num = 11
      )
    )
    expect_type(sc.out.cfg, "list")
    expect_true(length(sc.out.cfg) > 0, info = cfg$label)
  }

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
