library(npRmpi)

test_that("bounded cv.ls remains finite for gaussian order 2 and 4", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(), add = TRUE)

  set.seed(20260224)
  x <- runif(80)
  dat <- data.frame(x = x)

  bw2 <- npudensbw(
    dat = dat,
    bwmethod = "cv.ls",
    bwtype = "fixed",
    ckertype = "gaussian",
    ckerorder = 2L,
    ckerbound = "range",
    nmulti = 1
  )
  bw4 <- npudensbw(
    dat = dat,
    bwmethod = "cv.ls",
    bwtype = "fixed",
    ckertype = "gaussian",
    ckerorder = 4L,
    ckerbound = "range",
    nmulti = 1
  )

  expect_true(is.finite(as.numeric(bw2$bw[1])))
  expect_true(is.finite(as.numeric(bw2$fval)))
  expect_true(is.finite(as.numeric(bw4$bw[1])))
  expect_true(is.finite(as.numeric(bw4$fval)))
})

test_that("bounded unconditional cv.ls scalar quadrature supports generalized and adaptive NN bwtypes", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  set.seed(20260421)
  dat <- data.frame(x = runif(48L))

  bw.gnn <- npudensbw(
    dat = dat,
    bwmethod = "cv.ls",
    bwtype = "generalized_nn",
    ckertype = "epanechnikov",
    ckerorder = 6L,
    ckerbound = "range",
    nmulti = 1
  )

  bw.ad <- npudensbw(
    dat = dat,
    bwmethod = "cv.ls",
    bwtype = "adaptive_nn",
    ckertype = "gaussian",
    ckerorder = 8L,
    ckerbound = "range",
    nmulti = 1
  )

  expect_true(is.finite(as.numeric(bw.gnn$bw[1])))
  expect_true(is.finite(as.numeric(bw.gnn$fval)))
  expect_true(is.finite(as.numeric(bw.ad$bw[1])))
  expect_true(is.finite(as.numeric(bw.ad$fval)))
})

test_that("bounded unconditional cv.ls admits mixed and bivariate continuous bounded data", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  set.seed(20260421)
  mixed <- data.frame(x = runif(24L), g = factor(sample(c("a", "b"), 24L, replace = TRUE)))
  multi <- data.frame(x1 = runif(24L), x2 = runif(24L))
  combo <- data.frame(
    x1 = runif(24L),
    x2 = runif(24L),
    u = factor(sample(c("a", "b"), 24L, replace = TRUE)),
    o = ordered(sample(1:3, 24L, replace = TRUE))
  )

  bw.mixed <- npudensbw(
    dat = mixed,
    bwmethod = "cv.ls",
    bwtype = "generalized_nn",
    ckerbound = "range",
    nmulti = 1
  )

  bw.multi <- npudensbw(
    dat = multi,
    bwmethod = "cv.ls",
    bwtype = "adaptive_nn",
    ckerbound = "range",
    nmulti = 1
  )

  bw.combo <- npudensbw(
    dat = combo,
    bwmethod = "cv.ls",
    bwtype = "generalized_nn",
    ckerbound = "range",
    nmulti = 1
  )

  expect_true(all(is.finite(as.numeric(bw.mixed$bw))))
  expect_true(is.finite(as.numeric(bw.mixed$fval)))
  expect_true(all(is.finite(as.numeric(bw.multi$bw))))
  expect_true(is.finite(as.numeric(bw.multi$fval)))
  expect_true(all(is.finite(as.numeric(bw.combo$bw))))
  expect_true(is.finite(as.numeric(bw.combo$fval)))
})

test_that("bounded conditional cv.ls remains finite for gaussian order 2 and 4", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(), add = TRUE)

  set.seed(20260224)
  n <- 70
  x <- runif(n)
  y <- runif(n)

  xdat <- data.frame(x = x)
  ydat <- data.frame(y = y)

  bw2 <- npcdensbw(
    xdat = xdat,
    ydat = ydat,
    bwmethod = "cv.ls",
    bwtype = "fixed",
    cxkertype = "gaussian",
    cykertype = "gaussian",
    cxkerorder = 2L,
    cykerorder = 2L,
    cxkerbound = "range",
    cykerbound = "range",
    nmulti = 1
  )
  bw4 <- npcdensbw(
    xdat = xdat,
    ydat = ydat,
    bwmethod = "cv.ls",
    bwtype = "fixed",
    cxkertype = "gaussian",
    cykertype = "gaussian",
    cxkerorder = 4L,
    cykerorder = 4L,
    cxkerbound = "range",
    cykerbound = "range",
    nmulti = 1
  )

  expect_true(all(is.finite(as.numeric(bw2$xbw))))
  expect_true(all(is.finite(as.numeric(bw2$ybw))))
  expect_true(is.finite(as.numeric(bw2$fval)))
  expect_true(all(is.finite(as.numeric(bw4$xbw))))
  expect_true(all(is.finite(as.numeric(bw4$ybw))))
  expect_true(is.finite(as.numeric(bw4$fval)))
})

test_that("bounded conditional cv.ls scalar quadrature supports generalized and adaptive NN bwtypes", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  set.seed(20260421)
  n <- 48L
  xdat <- data.frame(x = runif(n))
  ydat <- data.frame(y = runif(n))

  bw.gnn <- npcdensbw(
    xdat = xdat,
    ydat = ydat,
    bwmethod = "cv.ls",
    bwtype = "generalized_nn",
    cxkertype = "epanechnikov",
    cykertype = "gaussian",
    cxkerorder = 4L,
    cykerorder = 8L,
    cxkerbound = "range",
    cykerbound = "range",
    nmulti = 1
  )

  bw.ad <- npcdensbw(
    xdat = xdat,
    ydat = ydat,
    bwmethod = "cv.ls",
    bwtype = "adaptive_nn",
    cxkertype = "gaussian",
    cykertype = "epanechnikov",
    cxkerorder = 2L,
    cykerorder = 6L,
    cxkerbound = "range",
    cykerbound = "range",
    nmulti = 1
  )

  expect_true(all(is.finite(as.numeric(bw.gnn$xbw))))
  expect_true(all(is.finite(as.numeric(bw.gnn$ybw))))
  expect_true(is.finite(as.numeric(bw.gnn$fval)))
  expect_true(all(is.finite(as.numeric(bw.ad$xbw))))
  expect_true(all(is.finite(as.numeric(bw.ad$ybw))))
  expect_true(is.finite(as.numeric(bw.ad$fval)))
})

test_that("bounded conditional cv.ls admits mixed and bivariate continuous bounded responses", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  set.seed(20260421)
  n <- 32L
  xdat <- data.frame(x = runif(n))
  ymixed <- data.frame(y = runif(n), g = factor(sample(c("a", "b"), n, replace = TRUE)))
  ymulti <- data.frame(y1 = runif(n), y2 = runif(n))
  ycombo <- data.frame(
    y1 = runif(n),
    y2 = runif(n),
    u = factor(sample(c("a", "b"), n, replace = TRUE)),
    o = ordered(sample(1:3, n, replace = TRUE))
  )

  bw.mixed <- npcdensbw(
    xdat = xdat,
    ydat = ymixed,
    bwmethod = "cv.ls",
    bwtype = "generalized_nn",
    cykerbound = "range",
    nmulti = 1
  )

  bw.multi <- npcdensbw(
    xdat = xdat,
    ydat = ymulti,
    bwmethod = "cv.ls",
    bwtype = "adaptive_nn",
    cykerbound = "range",
    nmulti = 1
  )

  bw.combo <- npcdensbw(
    xdat = xdat,
    ydat = ycombo,
    bwmethod = "cv.ls",
    bwtype = "generalized_nn",
    cykerbound = "range",
    nmulti = 1
  )

  expect_true(all(is.finite(as.numeric(bw.mixed$xbw))))
  expect_true(all(is.finite(as.numeric(bw.mixed$ybw))))
  expect_true(is.finite(as.numeric(bw.mixed$fval)))
  expect_true(all(is.finite(as.numeric(bw.multi$xbw))))
  expect_true(all(is.finite(as.numeric(bw.multi$ybw))))
  expect_true(is.finite(as.numeric(bw.multi$fval)))
  expect_true(all(is.finite(as.numeric(bw.combo$xbw))))
  expect_true(all(is.finite(as.numeric(bw.combo$ybw))))
  expect_true(is.finite(as.numeric(bw.combo$fval)))
})

test_that("bounded cv.ls still fails fast beyond two continuous bounded variables", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  set.seed(20260421)
  n <- 24L
  xdat <- data.frame(x = runif(n))
  y3 <- data.frame(y1 = runif(n), y2 = runif(n), y3 = runif(n))
  d3 <- data.frame(x1 = runif(n), x2 = runif(n), x3 = runif(n))

  expect_error(
    npcdensbw(
      xdat = xdat,
      ydat = y3,
      bwmethod = "cv.ls",
      bwtype = "fixed",
      cykerbound = "range",
      nmulti = 1
    ),
    "supports up to two continuous response variables"
  )

  expect_error(
    npudensbw(
      dat = d3,
      bwmethod = "cv.ls",
      bwtype = "fixed",
      ckerbound = "range",
      nmulti = 1
    ),
    "supports up to two continuous variables"
  )
})

test_that("bounded conditional distribution cv.ls remains finite for gaussian order 2 and 4", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(), add = TRUE)

  set.seed(20260312)
  n <- 70
  x <- runif(n)
  y <- runif(n)

  xdat <- data.frame(x = x)
  ydat <- data.frame(y = y)

  bw2 <- npcdistbw(
    xdat = xdat,
    ydat = ydat,
    bwmethod = "cv.ls",
    bwtype = "fixed",
    cxkertype = "gaussian",
    cykertype = "gaussian",
    cxkerorder = 2L,
    cykerorder = 2L,
    cxkerbound = "range",
    cykerbound = "range",
    nmulti = 1
  )
  bw4 <- npcdistbw(
    xdat = xdat,
    ydat = ydat,
    bwmethod = "cv.ls",
    bwtype = "fixed",
    cxkertype = "gaussian",
    cykertype = "gaussian",
    cxkerorder = 4L,
    cykerorder = 4L,
    cxkerbound = "range",
    cykerbound = "range",
    nmulti = 1
  )

  expect_true(all(is.finite(as.numeric(bw2$xbw))))
  expect_true(all(is.finite(as.numeric(bw2$ybw))))
  expect_true(is.finite(as.numeric(bw2$fval)))
  expect_true(all(is.finite(as.numeric(bw4$xbw))))
  expect_true(all(is.finite(as.numeric(bw4$ybw))))
  expect_true(is.finite(as.numeric(bw4$fval)))
})

test_that("bounded unconditional cv.ls remains finite after conditional bounded selectors", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(), add = TRUE)

  set.seed(20260316)
  n <- 48
  x <- runif(n)
  y <- runif(n)

  xdat <- data.frame(x = x)
  ydat <- data.frame(y = y)

  bw.cd <- npcdensbw(
    xdat = xdat,
    ydat = ydat,
    bwmethod = "cv.ls",
    bwtype = "fixed",
    cxkertype = "gaussian",
    cykertype = "gaussian",
    cxkerorder = 2L,
    cykerorder = 2L,
    cxkerbound = "range",
    cykerbound = "range",
    nmulti = 1
  )
  expect_true(all(is.finite(as.numeric(bw.cd$xbw))))
  expect_true(all(is.finite(as.numeric(bw.cd$ybw))))

  bw.ud <- npudensbw(
    dat = xdat,
    bwmethod = "cv.ls",
    bwtype = "fixed",
    ckertype = "gaussian",
    ckerorder = 4L,
    ckerbound = "range",
    nmulti = 1
  )

  expect_true(is.finite(as.numeric(bw.ud$bw[1])))
  expect_true(is.finite(as.numeric(bw.ud$fval)))
})
