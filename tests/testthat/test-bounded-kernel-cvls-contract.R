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
