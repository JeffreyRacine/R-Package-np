library(np)

test_that("bounded cv.ls remains finite for gaussian order 2 and 4", {
  set.seed(20260224)
  x <- runif(80)
  dat <- data.frame(x = x)

  for (ord in c(2L, 4L)) {
    bw <- npudensbw(
      dat = dat,
      bwmethod = "cv.ls",
      bwtype = "fixed",
      ckertype = "gaussian",
      ckerorder = ord,
      ckerbound = "range",
      nmulti = 1
    )

    expect_true(is.finite(as.numeric(bw$bw[1])))
    expect_true(is.finite(as.numeric(bw$fval)))
  }
})

test_that("bounded unconditional cv.ls scalar quadrature supports generalized and adaptive NN bwtypes", {
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

test_that("bounded unconditional cv.ls fails fast for unsupported data shapes", {
  set.seed(20260421)
  multi <- data.frame(x1 = runif(24L), x2 = runif(24L))
  mixed <- data.frame(x = runif(24L), g = factor(sample(c("a", "b"), 24L, replace = TRUE)))

  expect_error(
    npudensbw(
      dat = multi,
      bwmethod = "cv.ls",
      bwtype = "adaptive_nn",
      ckerbound = "range",
      nmulti = 1
    ),
    "bounded npudens cv\\.ls currently supports only one continuous variable"
  )

  expect_error(
    npudensbw(
      dat = mixed,
      bwmethod = "cv.ls",
      bwtype = "generalized_nn",
      ckerbound = "range",
      nmulti = 1
    ),
    "bounded npudens cv\\.ls currently supports only one continuous variable"
  )
})

test_that("bounded conditional cv.ls remains finite for gaussian order 2 and 4", {
  set.seed(20260224)
  n <- 70
  x <- runif(n)
  y <- runif(n)

  xdat <- data.frame(x = x)
  ydat <- data.frame(y = y)

  for (ord in c(2L, 4L)) {
    bw <- npcdensbw(
      xdat = xdat,
      ydat = ydat,
      bwmethod = "cv.ls",
      bwtype = "fixed",
      cxkertype = "gaussian",
      cykertype = "gaussian",
      cxkerorder = ord,
      cykerorder = ord,
      cxkerbound = "range",
      cykerbound = "range",
      nmulti = 1
    )

    expect_true(all(is.finite(as.numeric(bw$xbw))))
    expect_true(all(is.finite(as.numeric(bw$ybw))))
    expect_true(is.finite(as.numeric(bw$fval)))
  }
})

test_that("bounded conditional distribution cv.ls remains finite for gaussian order 2 and 4", {
  set.seed(20260224)
  n <- 70
  x <- runif(n)
  y <- runif(n)

  xdat <- data.frame(x = x)
  ydat <- data.frame(y = y)

  for (ord in c(2L, 4L)) {
    bw <- npcdistbw(
      xdat = xdat,
      ydat = ydat,
      bwmethod = "cv.ls",
      bwtype = "fixed",
      cxkertype = "gaussian",
      cykertype = "gaussian",
      cxkerorder = ord,
      cykerorder = ord,
      cxkerbound = "range",
      cykerbound = "range",
      nmulti = 1
    )

    expect_true(all(is.finite(as.numeric(bw$xbw))))
    expect_true(all(is.finite(as.numeric(bw$ybw))))
    expect_true(is.finite(as.numeric(bw$fval)))
  }
})

test_that("bounded unconditional cv.ls remains finite after conditional bounded selectors", {
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

test_that("bounded conditional cv.ls scalar quadrature supports generalized and adaptive NN bwtypes", {
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

test_that("bounded conditional cv.ls fails fast for unsupported response shapes", {
  set.seed(20260421)
  n <- 32L
  xdat <- data.frame(x = runif(n))
  ymulti <- data.frame(y1 = runif(n), y2 = runif(n))
  ymixed <- data.frame(y = runif(n), g = factor(sample(c("a", "b"), n, replace = TRUE)))

  expect_error(
    npcdensbw(
      xdat = xdat,
      ydat = ymulti,
      bwmethod = "cv.ls",
      bwtype = "adaptive_nn",
      cykerbound = "range",
      nmulti = 1
    ),
    "bounded response cv\\.ls currently supports only one continuous response variable"
  )

  expect_error(
    npcdensbw(
      xdat = xdat,
      ydat = ymixed,
      bwmethod = "cv.ls",
      bwtype = "generalized_nn",
      cykerbound = "range",
      nmulti = 1
    ),
    "bounded response cv\\.ls currently supports only one continuous response variable"
  )
})
