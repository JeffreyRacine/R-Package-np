test_that("semihat wrappers preserve infinite-bound parity and finite-bound evaluation checks", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(), add = TRUE)

  set.seed(20260320)
  n <- 90L
  x <- seq(0, 1, length.out = n)
  z <- seq(0, 1, length.out = n)
  y <- sin(2 * pi * z) + 0.8 * x + rnorm(n, sd = 0.03)

  tx <- data.frame(x = x)
  tz <- data.frame(z = z)
  ex <- data.frame(x = c(0, 0.5, 1))
  ez <- data.frame(z = c(0, 0.5, 1))
  ex.bad <- data.frame(x = c(-0.01, 0.5, 1.01))
  ez.bad <- data.frame(z = c(-0.01, 0.5, 1.01))

  sibw.none <- suppressWarnings(npindexbw(
    xdat = tx,
    ydat = y,
    bws = c(1, 0.4),
    bandwidth.compute = FALSE,
    regtype = "lc",
    bwtype = "fixed",
    ckertype = "gaussian",
    ckerorder = 2L,
    ckerbound = "none"
  ))
  sibw.inf <- suppressWarnings(npindexbw(
    xdat = tx,
    ydat = y,
    bws = c(1, 0.4),
    bandwidth.compute = FALSE,
    regtype = "lc",
    bwtype = "fixed",
    ckertype = "gaussian",
    ckerorder = 2L,
    ckerbound = "fixed",
    ckerlb = -Inf,
    ckerub = Inf
  ))
  sibw.fixed <- suppressWarnings(npindexbw(
    xdat = tx,
    ydat = y,
    bws = c(1, 0.4),
    bandwidth.compute = FALSE,
    regtype = "lc",
    bwtype = "fixed",
    ckertype = "gaussian",
    ckerorder = 2L,
    ckerbound = "fixed",
    ckerlb = 0,
    ckerub = 1
  ))

  si.none <- npindexhat(bws = sibw.none, txdat = tx, exdat = ex, y = y, s = 0L, output = "apply")
  si.inf <- npindexhat(bws = sibw.inf, txdat = tx, exdat = ex, y = y, s = 0L, output = "apply")
  si.fixed <- npindexhat(bws = sibw.fixed, txdat = tx, exdat = ex, y = y, s = 0L, output = "apply")

  expect_equal(as.vector(si.none), as.vector(si.inf), tolerance = 1e-10)
  expect_equal(as.vector(si.none), as.vector(si.fixed), tolerance = 1e-10)
  expect_error(
    npindexhat(bws = sibw.fixed, txdat = tx, exdat = ex.bad, y = y, s = 0L, output = "apply"),
    "Invalid bounds for 'ckerbound'|Evaluation data violate 'ckerbound' bounds"
  )

  scbw.none <- npscoefbw(
    xdat = tx,
    zdat = tz,
    ydat = y,
    bws = 0.4,
    bandwidth.compute = FALSE,
    regtype = "lc",
    bwtype = "fixed",
    ckertype = "gaussian",
    ckerorder = 2L,
    ckerbound = "none"
  )
  scbw.inf <- npscoefbw(
    xdat = tx,
    zdat = tz,
    ydat = y,
    bws = 0.4,
    bandwidth.compute = FALSE,
    regtype = "lc",
    bwtype = "fixed",
    ckertype = "gaussian",
    ckerorder = 2L,
    ckerbound = "fixed",
    ckerlb = -Inf,
    ckerub = Inf
  )
  scbw.fixed <- npscoefbw(
    xdat = tx,
    zdat = tz,
    ydat = y,
    bws = 0.4,
    bandwidth.compute = FALSE,
    regtype = "lc",
    bwtype = "fixed",
    ckertype = "gaussian",
    ckerorder = 2L,
    ckerbound = "fixed",
    ckerlb = 0,
    ckerub = 1
  )

  sc.none <- npscoefhat(
    bws = scbw.none,
    txdat = tx,
    tzdat = tz,
    exdat = ex,
    ezdat = ez,
    y = y,
    output = "apply",
    iterate = FALSE
  )
  sc.inf <- npscoefhat(
    bws = scbw.inf,
    txdat = tx,
    tzdat = tz,
    exdat = ex,
    ezdat = ez,
    y = y,
    output = "apply",
    iterate = FALSE
  )
  sc.fixed <- npscoefhat(
    bws = scbw.fixed,
    txdat = tx,
    tzdat = tz,
    exdat = ex,
    ezdat = ez,
    y = y,
    output = "apply",
    iterate = FALSE
  )

  expect_equal(as.vector(sc.none), as.vector(sc.inf), tolerance = 1e-10)
  expect_equal(as.vector(sc.none), as.vector(sc.fixed), tolerance = 1e-10)
  expect_error(
    npscoefhat(
      bws = scbw.fixed,
      txdat = tx,
      tzdat = tz,
      exdat = ex,
      ezdat = ez.bad,
      y = y,
      output = "apply",
      iterate = FALSE
    ),
    "Invalid bounds for 'ckerbound'|Evaluation data violate 'ckerbound' bounds"
  )

  plbw.none <- npplregbw(
    xdat = tx,
    zdat = tz,
    ydat = y,
    bws = matrix(c(0.4, 0.4), nrow = 2L, ncol = 1L),
    bandwidth.compute = FALSE,
    regtype = "lc",
    bwtype = "fixed",
    ckertype = "gaussian",
    ckerorder = 2L,
    ckerbound = "none"
  )
  plbw.inf <- npplregbw(
    xdat = tx,
    zdat = tz,
    ydat = y,
    bws = matrix(c(0.4, 0.4), nrow = 2L, ncol = 1L),
    bandwidth.compute = FALSE,
    regtype = "lc",
    bwtype = "fixed",
    ckertype = "gaussian",
    ckerorder = 2L,
    ckerbound = "fixed",
    ckerlb = -Inf,
    ckerub = Inf
  )
  plbw.fixed <- npplregbw(
    xdat = tx,
    zdat = tz,
    ydat = y,
    bws = matrix(c(0.4, 0.4), nrow = 2L, ncol = 1L),
    bandwidth.compute = FALSE,
    regtype = "lc",
    bwtype = "fixed",
    ckertype = "gaussian",
    ckerorder = 2L,
    ckerbound = "fixed",
    ckerlb = 0,
    ckerub = 1
  )

  pl.none <- npplreghat(
    bws = plbw.none,
    txdat = tx,
    tzdat = tz,
    exdat = ex,
    ezdat = ez,
    y = y,
    output = "apply"
  )
  pl.inf <- npplreghat(
    bws = plbw.inf,
    txdat = tx,
    tzdat = tz,
    exdat = ex,
    ezdat = ez,
    y = y,
    output = "apply"
  )
  pl.fixed <- npplreghat(
    bws = plbw.fixed,
    txdat = tx,
    tzdat = tz,
    exdat = ex,
    ezdat = ez,
    y = y,
    output = "apply"
  )

  expect_equal(as.vector(pl.none), as.vector(pl.inf), tolerance = 1e-10)
  expect_equal(as.vector(pl.none), as.vector(pl.fixed), tolerance = 1e-10)
  expect_error(
    npplreghat(
      bws = plbw.fixed,
      txdat = tx,
      tzdat = tz,
      exdat = ex,
      ezdat = ez.bad,
      y = y,
      output = "apply"
    ),
    "Invalid bounds for 'ckerbound'|Evaluation data violate 'ckerbound' bounds"
  )
})

test_that("bounded nonfixed semihat support is widened for pl, index, and scoef", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(), add = TRUE)

  set.seed(20260325)
  n <- 48L
  x1 <- runif(n)
  x2 <- runif(n)
  z <- runif(n)
  y_index <- sin(2 * pi * (0.7 * x1 - 0.4 * x2)) + rnorm(n, sd = 0.05)
  y_pl <- 0.8 * x1 + cos(2 * pi * z) + rnorm(n, sd = 0.05)
  y_sc <- (1 + sin(2 * pi * z)) * x1 + rnorm(n, sd = 0.05)

  tx_index <- data.frame(x1 = x1, x2 = x2)
  tx1 <- data.frame(x = x1)
  tz <- data.frame(z = z)

  bw.idx.gnn <- npindexbw(
    xdat = tx_index,
    ydat = y_index,
    method = "ichimura",
    bwtype = "generalized_nn",
    ckerbound = "range",
    nmulti = 1,
    itmax = 40,
    tol = 0.1
  )
  bw.idx.adapt <- npindexbw(
    xdat = tx_index,
    ydat = y_index,
    method = "ichimura",
    bwtype = "adaptive_nn",
    ckerbound = "range",
    nmulti = 1,
    itmax = 40,
    tol = 0.1
  )
  fit.idx.gnn <- npindex(bws = bw.idx.gnn, txdat = tx_index, tydat = y_index)
  fit.idx.adapt <- npindex(bws = bw.idx.adapt, txdat = tx_index, tydat = y_index)

  expect_true(all(is.finite(as.numeric(bw.idx.gnn$beta))))
  expect_true(is.finite(as.numeric(bw.idx.gnn$bw)))
  expect_true(all(is.finite(as.numeric(fit.idx.gnn$mean))))
  expect_true(all(is.finite(as.numeric(bw.idx.adapt$beta))))
  expect_true(is.finite(as.numeric(bw.idx.adapt$bw)))
  expect_true(all(is.finite(as.numeric(fit.idx.adapt$mean))))

  bw.pl.gnn <- npplregbw(
    xdat = tx1,
    ydat = y_pl,
    zdat = tz,
    regtype = "lc",
    bwtype = "generalized_nn",
    ckerbound = "range",
    nmulti = 1,
    itmax = 40,
    tol = 0.1
  )
  bw.pl.adapt <- npplregbw(
    xdat = tx1,
    ydat = y_pl,
    zdat = tz,
    regtype = "lc",
    bwtype = "adaptive_nn",
    ckerbound = "range",
    nmulti = 1,
    itmax = 40,
    tol = 0.1
  )
  fit.pl.gnn <- npplreg(bws = bw.pl.gnn, txdat = tx1, tydat = y_pl, tzdat = tz)
  fit.pl.adapt <- npplreg(bws = bw.pl.adapt, txdat = tx1, tydat = y_pl, tzdat = tz)

  expect_true(all(is.finite(as.numeric(unlist(lapply(bw.pl.gnn$bw, function(obj) obj$bw))))))
  expect_true(all(is.finite(as.numeric(unlist(lapply(bw.pl.adapt$bw, function(obj) obj$bw))))))
  expect_true(all(is.finite(as.numeric(fit.pl.gnn$mean))))
  expect_true(all(is.finite(as.numeric(fit.pl.adapt$mean))))

  bw.sc.gnn <- npscoefbw(
    xdat = tx1,
    ydat = y_sc,
    zdat = tz,
    regtype = "lc",
    bwtype = "generalized_nn",
    ckerbound = "range",
    nmulti = 1,
    itmax = 40,
    tol = 0.1
  )
  bw.sc.adapt <- npscoefbw(
    xdat = tx1,
    ydat = y_sc,
    zdat = tz,
    regtype = "lc",
    bwtype = "adaptive_nn",
    ckerbound = "range",
    nmulti = 1,
    itmax = 40,
    tol = 0.1
  )
  fit.sc.gnn <- npscoef(bws = bw.sc.gnn, txdat = tx1, tydat = y_sc, tzdat = tz, iterate = FALSE)
  fit.sc.adapt <- npscoef(bws = bw.sc.adapt, txdat = tx1, tydat = y_sc, tzdat = tz, iterate = FALSE)

  expect_true(all(is.finite(as.numeric(bw.sc.gnn$bw))))
  expect_true(all(is.finite(as.numeric(fit.sc.gnn$mean))))
  expect_true(all(is.finite(as.numeric(bw.sc.adapt$bw))))
  expect_true(all(is.finite(as.numeric(fit.sc.adapt$mean))))
})
