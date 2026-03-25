library(np)

test_that("fixed +/-Inf bounds are parity-equivalent to none at fixed bandwidth", {
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

test_that("bounded generalized_nn is available for certified core public routes", {
  set.seed(20260224)
  x <- runif(36)
  y <- cos(2 * pi * x) + rnorm(36, sd = 0.05)
  xy <- data.frame(x = x)
  yy <- data.frame(y = runif(36))

  bw.ud.cvml <- npudensbw(
    dat = xy,
    bwmethod = "cv.ml",
    bwtype = "generalized_nn",
    ckerbound = "range",
    nmulti = 1
  )
  bw.ud.cvls <- npudensbw(
    dat = xy,
    bwmethod = "cv.ls",
    bwtype = "generalized_nn",
    ckerbound = "range",
    nmulti = 1
  )
  fit.ud <- npudens(bws = bw.ud.cvls, tdat = xy)

  bw.dist <- npudistbw(
    dat = xy,
    bwmethod = "cv.cdf",
    bwtype = "generalized_nn",
    ckerbound = "range",
    nmulti = 1
  )
  fit.dist <- npudist(bws = bw.dist, tdat = xy)

  bw.reg <- npregbw(
    xdat = xy,
    ydat = y,
    regtype = "ll",
    bwmethod = "cv.ls",
    bwtype = "generalized_nn",
    ckerbound = "range",
    nmulti = 1
  )
  fit.reg <- npreg(bws = bw.reg, txdat = xy, tydat = y)

  bw.cd.cvml <- npcdensbw(
    xdat = xy,
    ydat = yy,
    bwmethod = "cv.ml",
    bwtype = "generalized_nn",
    cxkerbound = "range",
    cykerbound = "range",
    nmulti = 1
  )
  bw.cd.cvls <- npcdensbw(
    xdat = xy,
    ydat = yy,
    bwmethod = "cv.ls",
    bwtype = "generalized_nn",
    cxkerbound = "range",
    cykerbound = "range",
    nmulti = 1
  )
  fit.cd <- npcdens(bws = bw.cd.cvls, txdat = xy, tydat = yy)

  bw.cdist <- npcdistbw(
    xdat = xy,
    ydat = yy,
    bwmethod = "cv.ls",
    bwtype = "generalized_nn",
    cxkerbound = "range",
    cykerbound = "range",
    nmulti = 1
  )
  fit.cdist <- npcdist(bws = bw.cdist, txdat = xy, tydat = yy)

  expect_true(all(is.finite(as.numeric(bw.ud.cvml$bw))))
  expect_true(is.finite(bw.ud.cvml$fval))
  expect_true(all(is.finite(as.numeric(bw.ud.cvls$bw))))
  expect_true(is.finite(bw.ud.cvls$fval))
  expect_true(all(is.finite(as.numeric(fit.ud$dens))))

  expect_true(all(is.finite(as.numeric(bw.dist$bw))))
  expect_true(is.finite(bw.dist$fval))
  expect_true(all(is.finite(as.numeric(fit.dist$dist))))

  expect_true(all(is.finite(as.numeric(bw.reg$bw))))
  expect_true(is.finite(bw.reg$fval))
  expect_true(all(is.finite(as.numeric(fit.reg$mean))))

  expect_true(all(is.finite(as.numeric(bw.cd.cvml$xbw))))
  expect_true(all(is.finite(as.numeric(bw.cd.cvml$ybw))))
  expect_true(is.finite(bw.cd.cvml$fval))
  expect_true(all(is.finite(as.numeric(bw.cd.cvls$xbw))))
  expect_true(all(is.finite(as.numeric(bw.cd.cvls$ybw))))
  expect_true(is.finite(bw.cd.cvls$fval))
  expect_true(all(is.finite(as.numeric(fit.cd$condens))))

  expect_true(all(is.finite(as.numeric(bw.cdist$xbw))))
  expect_true(all(is.finite(as.numeric(bw.cdist$ybw))))
  expect_true(is.finite(bw.cdist$fval))
  expect_true(all(is.finite(as.numeric(fit.cdist$condist))))
})

test_that("npplreg, npindex, and npscoef bounded generalized_nn are available", {
  set.seed(20260224)
  x <- runif(32)
  y <- cos(2 * pi * x)
  y_sc <- (1 + sin(2 * pi * x)) * x + rnorm(32, sd = 0.05)
  xy <- data.frame(x = x)

  bw.pl <- npplregbw(
    xdat = xy,
    ydat = y,
    zdat = xy,
    bwmethod = "cv.ls",
    bwtype = "generalized_nn",
    ckerbound = "range",
    nmulti = 1
  )
  fit.pl <- npplreg(bws = bw.pl, txdat = xy, tydat = y, tzdat = xy)

  bw.index <- npindexbw(
    xdat = data.frame(x1 = x, x2 = x^2),
    ydat = y,
    method = "ichimura",
    bwtype = "generalized_nn",
    ckerbound = "range",
    nmulti = 1
  )
  fit.index <- npindex(bws = bw.index, txdat = data.frame(x1 = x, x2 = x^2), tydat = y)
  bw.sc <- npscoefbw(
    xdat = xy,
    ydat = y_sc,
    zdat = xy,
    bwtype = "generalized_nn",
    ckerbound = "range",
    nmulti = 1
  )
  fit.sc <- npscoef(bws = bw.sc, txdat = xy, tydat = y_sc, tzdat = xy, iterate = FALSE)

  expect_true(all(is.finite(as.numeric(unlist(lapply(bw.pl$bw, function(obj) obj$bw))))))
  expect_true(all(is.finite(as.numeric(fit.pl$mean))))
  expect_true(all(is.finite(as.numeric(bw.index$beta))))
  expect_true(is.finite(as.numeric(bw.index$bw)))
  expect_true(all(is.finite(as.numeric(fit.index$mean))))
  expect_true(all(is.finite(as.numeric(bw.sc$bw))))
  expect_true(all(is.finite(as.numeric(fit.sc$mean))))
})

test_that("evaluation support violations are caught before native execution", {
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
