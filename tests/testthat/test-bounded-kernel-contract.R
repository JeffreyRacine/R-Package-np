test_that("fixed +/-Inf bounds are parity-equivalent to none at fixed bandwidth", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  set.seed(20260224)
  x <- runif(80)
  dat <- data.frame(x = x)

  bw.none <- npudensbw(
    dat = dat,
    bws = 0.2,
    bandwidth.compute = FALSE,
    ckerbound = "none"
  )

  bw.inf <- npudensbw(
    dat = dat,
    bws = 0.2,
    bandwidth.compute = FALSE,
    ckerbound = "fixed",
    ckerlb = -Inf,
    ckerub = Inf
  )

  fit.none <- npudens(bws = bw.none, tdat = dat)
  fit.inf <- npudens(bws = bw.inf, tdat = dat)

  expect_equal(as.numeric(fit.none$dens), as.numeric(fit.inf$dens), tolerance = 1e-12)
})

test_that("scalar fixed bounds recycle over multiple continuous variables", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  set.seed(20260224)
  dat <- data.frame(x1 = runif(64), x2 = runif(64))

  bw.scalar <- npudensbw(
    dat = dat,
    bws = c(0.2, 0.25),
    bandwidth.compute = FALSE,
    ckerbound = "fixed",
    ckerlb = 0,
    ckerub = 1
  )

  bw.vector <- npudensbw(
    dat = dat,
    bws = c(0.2, 0.25),
    bandwidth.compute = FALSE,
    ckerbound = "fixed",
    ckerlb = c(0, 0),
    ckerub = c(1, 1)
  )

  fit.scalar <- npudens(bws = bw.scalar, tdat = dat)
  fit.vector <- npudens(bws = bw.vector, tdat = dat)

  expect_equal(as.numeric(fit.scalar$dens), as.numeric(fit.vector$dens), tolerance = 1e-12)
})

test_that("invalid fixed bounds are rejected with clear diagnostics", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  set.seed(20260224)
  x <- runif(50)
  dat <- data.frame(x = x)

  expect_error(
    npudensbw(
      dat = dat,
      bws = 0.2,
      bandwidth.compute = FALSE,
      ckerbound = "fixed",
      ckerlb = 1,
      ckerub = 1
    ),
    "lower < upper"
  )

  expect_error(
    npudensbw(
      dat = dat,
      bws = 0.2,
      bandwidth.compute = FALSE,
      ckerbound = "fixed",
      ckerlb = 0.2,
      ckerub = 0.8
    ),
    "Violations: x"
  )
})

test_that("finite bounds require fixed bwtype", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  set.seed(20260224)
  x <- runif(70)

  expect_error(
    npudensbw(
      dat = data.frame(x = x),
      bwmethod = "cv.ml",
      bwtype = "generalized_nn",
      ckerbound = "range",
      nmulti = 1
    ),
    "finite continuous kernel bounds require bwtype = \"fixed\""
  )
})

test_that("evaluation support violations are caught before native execution", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  set.seed(20260224)
  x <- runif(70)
  dat <- data.frame(x = x)

  bw <- npudensbw(
    dat = dat,
    bws = 0.2,
    bandwidth.compute = FALSE,
    ckerbound = "range"
  )

  expect_error(
    npudens(
      bws = bw,
      tdat = dat,
      edat = data.frame(x = max(x) + 0.05)
    ),
    "x >"
  )
})

test_that("predict paths enforce bounded eval checks with variable diagnostics", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  set.seed(20260224)
  x <- runif(80)
  y <- runif(80)
  dat.x <- data.frame(x = x)
  dat.y <- data.frame(y = y)

  bw.den <- npudensbw(
    dat = dat.x,
    bws = 0.2,
    bandwidth.compute = FALSE,
    ckerbound = "range"
  )
  fit.den <- npudens(bws = bw.den, tdat = dat.x)
  expect_error(
    predict(fit.den, edat = data.frame(x = max(x) + 0.01)),
    "Evaluation data violate 'ckerbound' bounds: x >"
  )

  bw.dist <- npudistbw(
    dat = dat.x,
    bws = 0.2,
    bandwidth.compute = FALSE,
    ckerbound = "range"
  )
  fit.dist <- npudist(bws = bw.dist, tdat = dat.x)
  expect_error(
    predict(fit.dist, edat = data.frame(x = max(x) + 0.01)),
    "Evaluation data violate 'ckerbound' bounds: x >"
  )

  bw.reg <- npregbw(
    xdat = dat.x,
    ydat = y,
    bws = 0.2,
    bandwidth.compute = FALSE,
    ckerbound = "range"
  )
  fit.reg <- npreg(bws = bw.reg, txdat = dat.x, tydat = y)
  expect_error(
    predict(fit.reg, exdat = data.frame(x = max(x) + 0.01)),
    "Evaluation data violate 'ckerbound' bounds: x >"
  )

  bw.cd <- npcdensbw(
    xdat = dat.x,
    ydat = dat.y,
    bws = c(0.2, 0.2),
    bandwidth.compute = FALSE,
    cxkerbound = "range",
    cykerbound = "range"
  )
  fit.cd <- npcdens(bws = bw.cd, txdat = dat.x, tydat = dat.y)
  expect_error(
    predict(fit.cd, exdat = data.frame(x = max(x) + 0.01), eydat = data.frame(y = mean(y))),
    "Evaluation data violate 'cxkerbound' bounds: x >"
  )
  expect_error(
    predict(fit.cd, exdat = data.frame(x = mean(x)), eydat = data.frame(y = max(y) + 0.01)),
    "Evaluation data violate 'cykerbound' bounds: y >"
  )

  bw.cdist <- npcdistbw(
    xdat = dat.x,
    ydat = dat.y,
    bws = c(0.2, 0.2),
    bandwidth.compute = FALSE,
    cxkerbound = "range",
    cykerbound = "range"
  )
  fit.cdist <- npcdist(bws = bw.cdist, txdat = dat.x, tydat = dat.y)
  expect_error(
    predict(fit.cdist, exdat = data.frame(x = max(x) + 0.01), eydat = data.frame(y = mean(y))),
    "Evaluation data violate 'cxkerbound' bounds: x >"
  )
  expect_error(
    predict(fit.cdist, exdat = data.frame(x = mean(x)), eydat = data.frame(y = max(y) + 0.01)),
    "Evaluation data violate 'cykerbound' bounds: y >"
  )
})
