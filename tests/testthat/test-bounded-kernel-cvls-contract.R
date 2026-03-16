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
