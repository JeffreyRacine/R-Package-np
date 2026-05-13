test_that("bandwidth predict methods forward object through bws", {
  old.options <- options(np.messages = FALSE)
  on.exit(options(old.options), add = TRUE)

  check_numeric_equal <- function(a, b, tolerance = 1e-10) {
    expect_equal(as.numeric(a), as.numeric(b), tolerance = tolerance)
  }

  set.seed(20260512)
  n <- 70
  x <- data.frame(x = runif(n))
  y <- sin(2 * pi * x$x) + rnorm(n, sd = 0.2)
  ey <- data.frame(y = y)
  z <- data.frame(z = runif(n))
  w <- data.frame(w = runif(n))
  nd <- data.frame(x = c(0.2, 0.5, 0.8))
  nd.yx <- data.frame(y = c(-0.4, 0, 0.4), x = c(0.2, 0.5, 0.8))
  nd.xw <- data.frame(x = c(0.2, 0.5, 0.8), w = c(0.1, 0.4, 0.7))

  bw.udens <- npudensbw(dat = x, bws = 0.25, bandwidth.compute = FALSE)
  check_numeric_equal(predict(bw.udens)$dens, npudens(bws = bw.udens)$dens)
  check_numeric_equal(
    predict(bw.udens, newdata = nd)$dens,
    npudens(bws = bw.udens, newdata = nd)$dens
  )

  bw.udist <- npudistbw(dat = x, bws = 0.25, bandwidth.compute = FALSE)
  check_numeric_equal(predict(bw.udist)$dist, npudist(bws = bw.udist)$dist)
  check_numeric_equal(
    predict(bw.udist, newdata = nd)$dist,
    npudist(bws = bw.udist, newdata = nd)$dist
  )

  bw.cdens <- npcdensbw(
    xdat = x,
    ydat = ey,
    bws = c(0.25, 0.2),
    bandwidth.compute = FALSE
  )
  check_numeric_equal(predict(bw.cdens)$condens, npcdens(bws = bw.cdens)$condens)
  check_numeric_equal(
    predict(bw.cdens, newdata = nd.yx)$condens,
    npcdens(bws = bw.cdens, newdata = nd.yx)$condens
  )

  bw.cdist <- npcdistbw(
    xdat = x,
    ydat = ey,
    bws = c(0.25, 0.2),
    bandwidth.compute = FALSE
  )
  check_numeric_equal(predict(bw.cdist)$condist, npcdist(bws = bw.cdist)$condist)
  check_numeric_equal(
    predict(bw.cdist, newdata = nd.yx)$condist,
    npcdist(bws = bw.cdist, newdata = nd.yx)$condist
  )

  bw.reg <- npregbw(xdat = x, ydat = y, bws = 0.25, bandwidth.compute = FALSE)
  check_numeric_equal(predict(bw.reg)$mean, npreg(bws = bw.reg)$mean)
  check_numeric_equal(
    predict(bw.reg, newdata = nd)$mean,
    npreg(bws = bw.reg, newdata = nd)$mean
  )

  bw.si <- npindexbw(
    xdat = data.frame(x = x$x, w = w$w),
    ydat = y,
    bws = c(0.25, 1, 0),
    bandwidth.compute = FALSE
  )
  check_numeric_equal(predict(bw.si)$mean, npindex(bws = bw.si)$mean)
  check_numeric_equal(
    predict(bw.si, newdata = nd.xw)$mean,
    npindex(bws = bw.si, newdata = nd.xw)$mean
  )

  semidat <- data.frame(
    y = y[seq_len(30)],
    x = x$x[seq_len(30)],
    z = z$z[seq_len(30)]
  )
  nd.xz <- data.frame(x = semidat$x[1:3], z = c(0.2, 0.5, 0.8))

  bw.pl <- npplregbw(y ~ x | z, data = semidat, nmulti = 1, tol = 0.1, ftol = 0.1)
  check_numeric_equal(predict(bw.pl)$mean, npplreg(bws = bw.pl)$mean)
  check_numeric_equal(
    predict(bw.pl, newdata = nd.xz)$mean,
    npplreg(bws = bw.pl, newdata = nd.xz)$mean
  )

  bw.sc <- npscoefbw(y ~ x | z, data = semidat, nmulti = 1, tol = 0.1, ftol = 0.1)
  check_numeric_equal(predict(bw.sc)$mean, npscoef(bws = bw.sc)$mean)
  check_numeric_equal(
    predict(bw.sc, newdata = nd.xz)$mean,
    npscoef(bws = bw.sc, newdata = nd.xz)$mean
  )
})
