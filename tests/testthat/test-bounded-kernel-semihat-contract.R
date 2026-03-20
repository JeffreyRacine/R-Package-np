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
