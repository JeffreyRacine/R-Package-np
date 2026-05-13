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

test_that("npindex and npindexhat preserve nearest-neighbor bwtype semantics", {
  set.seed(314162)
  n <- 80
  x1 <- runif(n)
  x2 <- runif(n)
  y <- sin(2 * (x1 + x2)) + rnorm(n, sd = 0.05)
  tx <- data.frame(x1 = x1, x2 = x2)
  ex <- tx[seq_len(25), , drop = FALSE]

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
      hy <- npindexhat(
        bws = bw,
        txdat = tx,
        exdat = ex,
        y = y,
        s = 0L,
        output = "apply"
      )

      expect_equal(as.vector(H %*% y), as.vector(fit$mean), tolerance = 1e-8, info = paste(bt, cfg$regtype))
      expect_equal(as.vector(hy), as.vector(fit$mean), tolerance = 1e-8, info = paste(bt, cfg$regtype))
    }
  }
})

test_that("npindexhat exact apply matches npindex on resampled nearest-neighbor lp fits", {
  set.seed(3292)
  n <- 45
  x1 <- runif(n)
  x2 <- runif(n)
  y <- sin(x1 + x2) + rnorm(n, sd = 0.08)
  tx <- data.frame(x1 = x1, x2 = x2)
  counts <- rmultinom(n = 3L, size = n, prob = rep.int(1 / n, n))

  for (bt in c("generalized_nn", "adaptive_nn")) {
    bw <- npindexbw(
      xdat = tx,
      ydat = y,
      bws = c(1, 1, 5L),
      bandwidth.compute = FALSE,
      regtype = "lp",
      basis = "tensor",
      degree = 2L,
      bwtype = bt
    )

    for (b in seq_len(ncol(counts))) {
      idx <- rep.int(seq_len(n), counts[, b])
      fit <- npindex(
        bws = bw,
        txdat = tx[idx, , drop = FALSE],
        tydat = y[idx],
        exdat = tx,
        gradients = FALSE
      )
      hy <- npindexhat(
        bws = bw,
        txdat = tx[idx, , drop = FALSE],
        exdat = tx,
        y = y[idx],
        output = "apply",
        s = 0L
      )

      expect_equal(as.vector(hy), as.vector(fit$mean), tolerance = 1e-8, info = paste(bt, "apply/public", b))
    }
  }
})

test_that("semihat apply mode matches core fits across bwtypes", {
  set.seed(20260308)
  n <- 70
  x <- runif(n)
  x2 <- runif(n)
  z <- runif(n)
  y <- (0.4 + x) * sin(2 * pi * z) + 0.3 * x2 + rnorm(n, sd = 0.04)

  tx.sc <- data.frame(x = x)
  tz.sc <- data.frame(z = z)
  ex.sc <- data.frame(x = seq(0.1, 0.9, length.out = 18))
  ez.sc <- data.frame(z = seq(0.1, 0.9, length.out = 18))

  tx.pl <- data.frame(x = x)
  tz.pl <- data.frame(z = z)
  ex.pl <- data.frame(x = seq(0.1, 0.9, length.out = 18))
  ez.pl <- data.frame(z = seq(0.1, 0.9, length.out = 18))

  tx.si <- data.frame(x1 = x, x2 = x2)
  ex.si <- tx.si[seq_len(18), , drop = FALSE]

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
      pl.args <- list(
        xdat = tx.pl,
        zdat = tz.pl,
        ydat = y,
        regtype = regtype,
        bwtype = bwtype,
        bandwidth.compute = FALSE,
        bws = matrix(rep(if (identical(bwtype, "fixed")) 0.18 else 9, 2L), nrow = 2L, ncol = 1L)
      )
      si.args <- list(
        xdat = tx.si,
        ydat = y,
        regtype = regtype,
        bwtype = bwtype,
        bandwidth.compute = FALSE,
        bws = c(1, 1, if (identical(bwtype, "fixed")) 0.18 else 9)
      )
      if (identical(regtype, "lp")) {
        sc.args$degree <- 1L
        sc.args$basis <- "glp"
        sc.args$bernstein.basis <- FALSE
        pl.args$degree <- 1L
        pl.args$basis <- "glp"
        pl.args$bernstein.basis <- FALSE
        si.args$degree <- 1L
        si.args$basis <- "glp"
        si.args$bernstein.basis <- FALSE
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
      sc.constraint <- npscoefhat(
        bws = sc.bw,
        txdat = tx.sc,
        tzdat = tz.sc,
        exdat = ex.sc,
        ezdat = ez.sc,
        y = y,
        output = "constraint",
        iterate = FALSE
      )
      expect_equal(
        as.vector(sc.apply),
        as.vector(sc.H %*% y),
        tolerance = 1e-10,
        info = paste("npscoefhat helper parity", regtype, bwtype)
      )
      expect_equal(
        sc.constraint,
        t(sc.H) * y,
        tolerance = 0,
        ignore_attr = TRUE,
        info = paste("npscoefhat constraint", regtype, bwtype)
      )
      if (identical(bwtype, "fixed")) {
        expect_equal(
          as.vector(sc.apply),
          as.vector(sc.fit$mean),
          tolerance = 1e-8,
          info = paste("npscoefhat core parity", regtype, bwtype)
        )
      }

      pl.bw <- do.call(npplregbw, pl.args)
      pl.fit <- npplreg(
        bws = pl.bw,
        txdat = tx.pl,
        tzdat = tz.pl,
        tydat = y,
        exdat = ex.pl,
        ezdat = ez.pl
      )
      pl.apply <- npplreghat(
        bws = pl.bw,
        txdat = tx.pl,
        tzdat = tz.pl,
        exdat = ex.pl,
        ezdat = ez.pl,
        y = y,
        output = "apply"
      )
      pl.H <- npplreghat(
        bws = pl.bw,
        txdat = tx.pl,
        tzdat = tz.pl,
        exdat = ex.pl,
        ezdat = ez.pl,
        output = "matrix"
      )
      pl.constraint <- npplreghat(
        bws = pl.bw,
        txdat = tx.pl,
        tzdat = tz.pl,
        exdat = ex.pl,
        ezdat = ez.pl,
        y = y,
        output = "constraint"
      )
      expect_equal(
        as.vector(pl.apply),
        as.vector(pl.fit$mean),
        tolerance = 1e-8,
        info = paste("npplreghat", regtype, bwtype)
      )
      expect_equal(
        pl.constraint,
        t(pl.H) * y,
        tolerance = 0,
        ignore_attr = TRUE,
        info = paste("npplreghat constraint", regtype, bwtype)
      )

      si.bw <- do.call(npindexbw, si.args)
      si.fit.mean <- npindex(
        bws = si.bw,
        txdat = tx.si,
        tydat = y,
        exdat = ex.si,
        gradients = FALSE
      )
      si.fit.grad <- npindex(
        bws = si.bw,
        txdat = tx.si,
        tydat = y,
        exdat = ex.si,
        gradients = TRUE
      )
      si.apply.mean <- npindexhat(
        bws = si.bw,
        txdat = tx.si,
        exdat = ex.si,
        y = y,
        output = "apply",
        s = 0L
      )
      si.H.mean <- npindexhat(
        bws = si.bw,
        txdat = tx.si,
        exdat = ex.si,
        output = "matrix",
        s = 0L
      )
      si.constraint.mean <- npindexhat(
        bws = si.bw,
        txdat = tx.si,
        exdat = ex.si,
        y = y,
        output = "constraint",
        s = 0L
      )
      si.apply.grad <- npindexhat(
        bws = si.bw,
        txdat = tx.si,
        exdat = ex.si,
        y = y,
        output = "apply",
        s = 1L
      )
      si.H.grad <- npindexhat(
        bws = si.bw,
        txdat = tx.si,
        exdat = ex.si,
        output = "matrix",
        s = 1L
      )
      si.constraint.grad <- npindexhat(
        bws = si.bw,
        txdat = tx.si,
        exdat = ex.si,
        y = y,
        output = "constraint",
        s = 1L
      )
      expect_equal(
        as.vector(si.apply.mean),
        as.vector(si.fit.mean$mean),
        tolerance = 1e-8,
        info = paste("npindexhat mean", regtype, bwtype)
      )
      expect_equal(
        as.vector(si.H.mean %*% y),
        as.vector(si.fit.mean$mean),
        tolerance = 1e-8,
        info = paste("npindexhat mean matrix", regtype, bwtype)
      )
      expect_equal(
        as.vector(si.H.mean %*% y),
        as.vector(si.apply.mean),
        tolerance = 1e-10,
        info = paste("npindexhat mean matrix/apply", regtype, bwtype)
      )
      expect_equal(
        si.constraint.mean,
        t(si.H.mean) * y,
        tolerance = 0,
        ignore_attr = TRUE,
        info = paste("npindexhat mean constraint", regtype, bwtype)
      )
      expect_equal(
        as.vector(si.apply.grad),
        as.vector(si.fit.grad$grad[, 1]),
        tolerance = 1e-8,
        info = paste("npindexhat grad", regtype, bwtype)
      )
      expect_equal(
        as.vector(si.H.grad %*% y),
        as.vector(si.fit.grad$grad[, 1]),
        tolerance = 1e-8,
        info = paste("npindexhat grad matrix", regtype, bwtype)
      )
      expect_equal(
        as.vector(si.H.grad %*% y),
        as.vector(si.apply.grad),
        tolerance = 1e-10,
        info = paste("npindexhat grad matrix/apply", regtype, bwtype)
      )
      expect_equal(
        si.constraint.grad,
        t(si.H.grad) * y,
        tolerance = 0,
        ignore_attr = TRUE,
        info = paste("npindexhat grad constraint", regtype, bwtype)
      )
    }
  }
})

test_that("npscoefbw adaptive-nn manual bandwidth must be integer support", {
  set.seed(20260309)
  n <- 40L
  x <- runif(n)
  z <- runif(n)
  y <- (0.4 + x) * sin(2 * pi * z) + rnorm(n, sd = 0.04)
  tx <- data.frame(x = x)
  tz <- data.frame(z = z)

  expect_error(
    npscoefbw(
      xdat = tx,
      zdat = tz,
      ydat = y,
      regtype = "lc",
      bwtype = "adaptive_nn",
      bandwidth.compute = FALSE,
      bws = 0.13
    ),
    "nearest-neighbor bandwidth must be an integer"
  )
})

test_that("npscoefhat selected adaptive-nn owner preserves integer support and helper parity", {
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
  round_half_to_even <- getFromNamespace(".np_round_half_to_even", "np")

  for (regtype in c("lc", "ll", "lp")) {
    bw_args <- list(
      xdat = tx,
      zdat = tz,
      ydat = y,
      regtype = regtype,
      bwtype = "adaptive_nn",
      nmulti = 1L
    )
    if (identical(regtype, "lp")) {
      bw_args$degree <- 1L
      bw_args$basis <- "glp"
      bw_args$bernstein.basis <- FALSE
    }

    bw <- do.call(npscoefbw, bw_args)
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

test_that("npindexbw lc selection no longer collapses nearest-neighbor bwtypes to fixed", {
  set.seed(314163)
  n <- 70
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  y <- x1 - x2 + rnorm(n, sd = 0.2)
  tx <- data.frame(x1 = x1, x2 = x2)

  bw.fixed <- npindexbw(xdat = tx, ydat = y, regtype = "lc", bwtype = "fixed", nmulti = 1)
  bw.gen <- npindexbw(xdat = tx, ydat = y, regtype = "lc", bwtype = "generalized_nn", nmulti = 1)
  bw.adap <- npindexbw(xdat = tx, ydat = y, regtype = "lc", bwtype = "adaptive_nn", nmulti = 1)

  fit.fixed <- npindex(bws = bw.fixed, txdat = tx, tydat = y)$mean
  fit.gen <- npindex(bws = bw.gen, txdat = tx, tydat = y)$mean
  fit.adap <- npindex(bws = bw.adap, txdat = tx, tydat = y)$mean

  expect_gt(max(abs(fit.fixed - fit.gen)), 1e-6)
  expect_gt(max(abs(fit.fixed - fit.adap)), 1e-6)
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
    "Invalid bounds for 'ckerbound'|Evaluation data violate 'ckerbound' bounds"
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
  expect_error(
    npindexhat(bws = structure(list(), class = "sibandwidth"),
               txdat = data.frame(x = x), foo = TRUE),
    "unused argument in npindexhat: 'foo'"
  )
  expect_error(
    npplreghat(bws = structure(list(), class = "plbandwidth"),
               txdat = data.frame(x = x), tzdat = data.frame(z = z), foo = TRUE),
    "unused argument in npplreghat: 'foo'"
  )
  expect_error(
    npscoefhat(bws = structure(list(), class = "scbandwidth"),
               txdat = data.frame(x = x), tzdat = data.frame(z = z), foo = TRUE),
    "unused argument in npscoefhat: 'foo'"
  )

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
      output = "data",
      errors = "bootstrap",
      bootstrap = "wild",
      B = 19
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
      output = "data",
      errors = "bootstrap",
      bootstrap = "wild",
      boot_control = np_boot_control(wild = "rademacher"),
      band = "pointwise",
      B = 19
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
        output = "data",
        errors = "bootstrap",
        bootstrap = "wild",
        B = 11
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
        output = "data",
        errors = "bootstrap",
        bootstrap = "wild",
        B = 9
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
      output = "data",
      errors = "bootstrap",
      bootstrap = "wild",
      B = 19
    )
  )
  expect_type(si.out, "list")
  expect_true(length(si.out) > 0)
  expect_true(is.matrix(si.out[[1]]$merr))
  expect_equal(ncol(si.out[[1]]$merr), 2L)
  expect_false(all(is.na(si.out[[1]]$merr)))
})
