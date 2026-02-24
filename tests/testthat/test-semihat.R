test_that("npscoefhat reproduces npscoef fitted values and supports matrix RHS", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  old.auto <- getOption("npRmpi.autodispatch", FALSE)
  on.exit(options(npRmpi.autodispatch = old.auto), add = TRUE)
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
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
})

test_that("npplreghat reproduces npplreg fitted values and supports matrix RHS", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  old.auto <- getOption("npRmpi.autodispatch", FALSE)
  on.exit(options(npRmpi.autodispatch = old.auto), add = TRUE)
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
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

test_that("npplreghat supports ll/lp with lp basis variants", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  old.auto <- getOption("npRmpi.autodispatch", FALSE)
  on.exit(options(npRmpi.autodispatch = old.auto), add = TRUE)
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
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
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
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
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
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


test_that("semihat validates class and scalar controls", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  old.auto <- getOption("npRmpi.autodispatch", FALSE)
  on.exit(options(npRmpi.autodispatch = old.auto), add = TRUE)
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
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

  scbw <- npscoefbw(
    xdat = x,
    zdat = z,
    ydat = y,
    bws = 0.2,
    bandwidth.compute = FALSE
  )
  expect_error(
    npscoefhat(
      bws = scbw,
      txdat = data.frame(x = x),
      tzdat = data.frame(z = z),
      ridge = 0
    ),
    "positive finite scalar"
  )
})

test_that("plot bootstrap supports wild for sc/pl/si bandwidth objects", {
  skip_on_cran()
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  old.auto <- getOption("npRmpi.autodispatch", FALSE)
  on.exit(options(npRmpi.autodispatch = old.auto), add = TRUE)
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
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
