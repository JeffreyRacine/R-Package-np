library(np)

test_that("predict aliases newdata to native eval args for default npreg/npudens/npudist/npindex", {
  set.seed(20260224)
  x <- runif(70)
  y <- rnorm(70)
  nd <- data.frame(x = c(0.1, 0.3, 0.7))

  bw.reg <- npregbw(
    xdat = data.frame(x = x),
    ydat = y,
    bws = 0.25,
    bandwidth.compute = FALSE
  )
  fit.reg <- npreg(bws = bw.reg, txdat = data.frame(x = x), tydat = y)
  expect_equal(
    as.numeric(predict(fit.reg, newdata = nd)),
    as.numeric(predict(fit.reg, exdat = nd)),
    tolerance = 1e-12
  )

  bw.den <- npudensbw(
    dat = data.frame(x = x),
    bws = 0.25,
    bandwidth.compute = FALSE
  )
  fit.den <- npudens(bws = bw.den, tdat = data.frame(x = x))
  expect_equal(
    as.numeric(predict(fit.den, newdata = nd)),
    as.numeric(predict(fit.den, edat = nd)),
    tolerance = 1e-12
  )

  bw.dist <- npudistbw(
    dat = data.frame(x = x),
    bws = 0.25,
    bandwidth.compute = FALSE
  )
  fit.dist <- npudist(bws = bw.dist, tdat = data.frame(x = x))
  expect_equal(
    as.numeric(predict(fit.dist, newdata = nd)),
    as.numeric(predict(fit.dist, edat = nd)),
    tolerance = 1e-12
  )

  x2 <- runif(70)
  nd.si <- data.frame(x = c(0.15, 0.35), x2 = c(0.4, 0.8))
  bw.si <- npindexbw(
    xdat = data.frame(x = x, x2 = x2),
    ydat = y,
    bws = c(0.25, 0.25, 1),
    bandwidth.compute = FALSE
  )
  fit.si <- npindex(
    bws = bw.si,
    txdat = data.frame(x = x, x2 = x2),
    tydat = y
  )
  expect_equal(
    as.numeric(predict(fit.si, newdata = nd.si)),
    as.numeric(predict(fit.si, exdat = nd.si)),
    tolerance = 1e-12
  )
})

test_that("predict aliases newdata to exdat/eydat for default npcdens/npcdist", {
  set.seed(20260224)
  x <- runif(60)
  y <- runif(60)
  nd <- data.frame(y = c(0.2, 0.5), x = c(0.1, 0.8))

  bw.cd <- npcdensbw(
    xdat = data.frame(x = x),
    ydat = data.frame(y = y),
    bws = c(0.25, 0.2),
    bandwidth.compute = FALSE
  )
  fit.cd <- npcdens(
    bws = bw.cd,
    txdat = data.frame(x = x),
    tydat = data.frame(y = y)
  )
  expect_equal(
    as.numeric(predict(fit.cd, newdata = nd)),
    as.numeric(predict(fit.cd, exdat = nd["x"], eydat = nd["y"])),
    tolerance = 1e-12
  )
  expect_error(
    predict(fit.cd, newdata = data.frame(x = c(0.1, 0.2))),
    "must include columns"
  )

  bw.cdist <- npcdistbw(
    xdat = data.frame(x = x),
    ydat = data.frame(y = y),
    bws = c(0.25, 0.2),
    bandwidth.compute = FALSE
  )
  fit.cdist <- npcdist(
    bws = bw.cdist,
    txdat = data.frame(x = x),
    tydat = data.frame(y = y)
  )
  expect_equal(
    as.numeric(predict(fit.cdist, newdata = nd)),
    as.numeric(predict(fit.cdist, exdat = nd["x"], eydat = nd["y"])),
    tolerance = 1e-12
  )
})

test_that("predict aliases newdata to exdat for default npqreg", {
  set.seed(20260512)
  x <- data.frame(x = runif(50))
  y <- rnorm(50)
  nd <- data.frame(x = c(0.2, 0.5, 0.8))

  bw <- npcdistbw(
    xdat = x,
    ydat = data.frame(y = y),
    bws = c(0.25, 0.25),
    bandwidth.compute = FALSE
  )
  fit <- npqreg(bws = bw, txdat = x, tydat = y, tau = 0.5)

  expect_equal(
    as.numeric(predict(fit, newdata = nd)),
    as.numeric(predict(fit, exdat = nd)),
    tolerance = 1e-12
  )

  nd.native <- data.frame(x = c(0.15, 0.35, 0.65))
  nd.standard <- data.frame(x = c(0.25, 0.55, 0.85))
  expect_equal(
    as.numeric(predict(fit, newdata = nd.standard, exdat = nd.native)),
    as.numeric(predict(fit, exdat = nd.native)),
    tolerance = 1e-12
  )
  expect_error(
    predict(fit, newdata = data.frame(z = c(0.2, 0.5, 0.8))),
    "newdata must contain columns"
  )
})

test_that("predict supports formula newdata and vector tau for npqreg", {
  set.seed(20260512)
  dat <- data.frame(
    y = sin(runif(70)) + rnorm(70, sd = 0.1),
    x = runif(70)
  )
  nd <- data.frame(x = c(0.2, 0.5, 0.8))

  bw <- npcdistbw(
    y ~ x,
    data = dat,
    bws = c(0.25, 0.25),
    bandwidth.compute = FALSE
  )
  fit <- npqreg(bws = bw, tau = c(0.25, 0.5), gradients = TRUE)

  pred <- predict(fit, newdata = nd)
  expect_equal(dim(pred), c(nrow(nd), 2L))

  se.pred <- predict(fit, newdata = nd, se.fit = TRUE)
  expect_named(se.pred, c("fit", "se.fit", "df", "residual.scale"))
  expect_equal(dim(se.pred$fit), c(nrow(nd), 2L))
  expect_equal(dim(se.pred$se.fit), c(nrow(nd), 2L))

  nd.native <- data.frame(x = c(0.15, 0.35))
  expect_equal(
    as.numeric(predict(fit, exdat = nd.native, newdata = nd)),
    as.numeric(predict(fit, newdata = nd.native)),
    tolerance = 1e-12
  )
  expect_error(
    predict(fit, se.fit = "yes"),
    "'se.fit' must be TRUE or FALSE",
    fixed = TRUE
  )
  x <- rep(0.1, nrow(nd))
  expect_error(
    npqreg(bws = bw, newdata = data.frame(z = x), tau = 0.5),
    "newdata must contain columns"
  )
  expect_error(
    predict(fit, newdata = data.frame(z = x)),
    "newdata must contain columns"
  )
})

test_that("predict supports conmode formula newdata without changing extraction behavior", {
  set.seed(20260512)
  dat <- data.frame(
    y = factor(rbinom(80, 1, 0.45)),
    x = runif(80),
    z = factor(rbinom(80, 1, 0.5))
  )
  bw <- npcdensbw(
    y ~ x + z,
    data = dat,
    bws = c(0.25, 0.3, 0.4),
    bandwidth.compute = FALSE
  )
  fit <- npconmode(bws = bw, probabilities = TRUE)
  fit.no.prob <- npconmode(bws = bw, probabilities = FALSE)
  nd <- dat[1:5, c("x", "z")]
  nd.with.y <- dat[1:5, c("y", "x", "z")]

  expect_equal(predict(fit), fit$conmode)
  expect_equal(predict(fit, type = "prob"), fit$probabilities)
  expect_error(
    predict(fit.no.prob, type = "prob"),
    "class probabilities are not stored"
  )

  ref <- npconmode(bws = fit$bws, newdata = nd)
  expect_equal(as.character(predict(fit, newdata = nd)), as.character(ref$conmode))

  ref.prob <- npconmode(bws = fit$bws, newdata = nd, probabilities = TRUE)
  expect_equal(
    unname(predict(fit, newdata = nd, type = "prob")),
    unname(ref.prob$probabilities),
    tolerance = 1e-12
  )

  ref.with.y <- npconmode(bws = fit$bws, newdata = nd.with.y)
  expect_equal(
    as.character(predict(fit, newdata = nd.with.y)),
    as.character(ref.with.y$conmode)
  )
  expect_error(
    predict(fit, newdata = nd, type = "prob", probabilities = FALSE),
    "requires probabilities=TRUE"
  )
  x <- rep(0.1, nrow(nd))
  expect_error(
    npconmode(bws = bw, newdata = data.frame(wrong = x)),
    "newdata must contain columns"
  )
  expect_error(
    predict(fit, newdata = data.frame(wrong = x)),
    "newdata must contain columns"
  )
})

test_that("predict supports conmode native evaluation and native precedence", {
  set.seed(20260512)
  xdat <- data.frame(
    x = runif(70),
    z = factor(rbinom(70, 1, 0.5))
  )
  ydat <- factor(rbinom(70, 1, 0.4))
  bw <- npcdensbw(
    xdat = xdat,
    ydat = data.frame(y = ydat),
    bws = c(0.25, 0.3, 0.4),
    bandwidth.compute = FALSE
  )
  fit <- npconmode(bws = bw, txdat = xdat, tydat = ydat, probabilities = TRUE)
  nd <- xdat[1:5, , drop = FALSE]
  nd.native <- xdat[6:10, , drop = FALSE]

  expect_equal(
    as.character(predict(fit, newdata = nd)),
    as.character(predict(fit, exdat = nd))
  )
  expect_equal(
    unname(predict(fit, newdata = nd, type = "prob")),
    unname(predict(fit, exdat = nd, type = "prob")),
    tolerance = 1e-12
  )
  expect_equal(
    as.character(predict(fit, newdata = nd, exdat = NULL)),
    as.character(predict(fit, newdata = nd))
  )
  expect_error(
    predict(fit, eydat = data.frame(y = ydat[1:5])),
    "eydat.*exdat"
  )
  expect_equal(
    as.character(predict(fit, newdata = nd, exdat = nd.native)),
    as.character(predict(fit, exdat = nd.native))
  )
})

test_that("predict aliases newdata to the explicit-evaluation slice route for npcdens", {
  set.seed(20260322)
  x <- runif(70, -1, 1)
  y <- sin(2 * pi * x) + rnorm(70, sd = 0.2)
  nd <- rbind(
    data.frame(y = c(-0.7, -0.15, 0.45), x = rep(-0.35, 3L)),
    data.frame(y = c(-0.25, 0.2, 0.75, 1.0), x = rep(0.4, 4L))
  )
  ctrl <- list(mode = "slice", slice.grid.size = 21L, slice.extend.factor = 0)

  bw.cd <- npcdensbw(
    xdat = data.frame(x = x),
    ydat = data.frame(y = y),
    bws = c(0.27, 0.21),
    bandwidth.compute = FALSE,
    regtype = "lp",
    degree = 3L
  )
  fit.cd <- npcdens(
    bws = bw.cd,
    txdat = data.frame(x = x),
    tydat = data.frame(y = y),
    proper = TRUE
  )

  expect_equal(
    as.numeric(predict(fit.cd, newdata = nd, proper.control = ctrl)),
    as.numeric(predict(
      fit.cd,
      exdat = nd["x"],
      eydat = nd["y"],
      proper.control = ctrl
    )),
    tolerance = 1e-10
  )
})

test_that("predict newdata is unchanged when proper apply='fitted' targets only fitted values", {
  set.seed(20260322)
  x <- runif(70, -1, 1)
  y <- sin(2 * pi * x) + rnorm(70, sd = 0.2)
  nd <- rbind(
    data.frame(y = c(-0.7, -0.15, 0.45), x = rep(-0.35, 3L)),
    data.frame(y = c(-0.25, 0.2, 0.75, 1.0), x = rep(0.4, 4L))
  )
  ctrl <- list(mode = "slice", apply = "fitted", slice.grid.size = 21L, slice.extend.factor = 0)

  bw.cd <- npcdensbw(
    xdat = data.frame(x = x),
    ydat = data.frame(y = y),
    bws = c(0.27, 0.21),
    bandwidth.compute = FALSE,
    regtype = "lp",
    degree = 3L
  )
  fit.cd <- npcdens(
    bws = bw.cd,
    txdat = data.frame(x = x),
    tydat = data.frame(y = y),
    proper = TRUE
  )

  expect_equal(
    as.numeric(predict(fit.cd, newdata = nd, proper = FALSE)),
    as.numeric(predict(fit.cd, newdata = nd, proper.control = ctrl)),
    tolerance = 1e-12
  )
})

test_that("predict aliases newdata to the explicit-evaluation slice route for npcdist", {
  set.seed(20260322)
  x <- runif(70, -1, 1)
  y <- sin(2 * pi * x) + rnorm(70, sd = 0.2)
  nd <- rbind(
    data.frame(y = c(-0.7, -0.15, 0.45), x = rep(-0.35, 3L)),
    data.frame(y = c(-0.25, 0.2, 0.75, 1.0), x = rep(0.4, 4L))
  )
  ctrl <- list(mode = "slice", slice.grid.size = 21L, slice.extend.factor = 0)

  bw.cd <- npcdistbw(
    xdat = data.frame(x = x),
    ydat = data.frame(y = y),
    bws = c(0.27, 0.21),
    bandwidth.compute = FALSE,
    regtype = "lp",
    degree = 3L
  )
  fit.cd <- npcdist(
    bws = bw.cd,
    txdat = data.frame(x = x),
    tydat = data.frame(y = y),
    proper = TRUE
  )

  expect_equal(
    as.numeric(predict(fit.cd, newdata = nd, proper.control = ctrl)),
    as.numeric(predict(
      fit.cd,
      exdat = nd["x"],
      eydat = nd["y"],
      proper.control = ctrl
    )),
    tolerance = 1e-10
  )
})

test_that("predict newdata for npcdist is unchanged when proper apply='fitted' targets only fitted values", {
  set.seed(20260322)
  x <- runif(70, -1, 1)
  y <- sin(2 * pi * x) + rnorm(70, sd = 0.2)
  nd <- rbind(
    data.frame(y = c(-0.7, -0.15, 0.45), x = rep(-0.35, 3L)),
    data.frame(y = c(-0.25, 0.2, 0.75, 1.0), x = rep(0.4, 4L))
  )
  ctrl <- list(mode = "slice", apply = "fitted", slice.grid.size = 21L, slice.extend.factor = 0)

  bw.cd <- npcdistbw(
    xdat = data.frame(x = x),
    ydat = data.frame(y = y),
    bws = c(0.27, 0.21),
    bandwidth.compute = FALSE,
    regtype = "lp",
    degree = 3L
  )
  fit.cd <- npcdist(
    bws = bw.cd,
    txdat = data.frame(x = x),
    tydat = data.frame(y = y),
    proper = TRUE
  )

  expect_equal(
    as.numeric(predict(fit.cd, newdata = nd, proper = FALSE)),
    as.numeric(predict(fit.cd, newdata = nd, proper.control = ctrl)),
    tolerance = 1e-12
  )
})

test_that("predict aliases newdata to exdat/ezdat for default npscoef/npplreg", {
  set.seed(20260227)
  n <- 80
  x <- runif(n)
  z <- runif(n)
  y <- sin(2 * pi * x) + 1.5 * z + rnorm(n, sd = 0.05)
  nd <- data.frame(x = c(0.2, 0.5, 0.8), z = c(0.1, 0.4, 0.9))

  bw.sc <- npscoefbw(
    xdat = data.frame(x = x),
    ydat = y,
    zdat = data.frame(z = z),
    bws = 0.25,
    bandwidth.compute = FALSE
  )
  fit.sc <- npscoef(
    bws = bw.sc,
    txdat = data.frame(x = x),
    tydat = y,
    tzdat = data.frame(z = z)
  )
  expect_equal(
    as.numeric(predict(fit.sc, newdata = nd)),
    as.numeric(predict(fit.sc, exdat = nd["x"], ezdat = nd["z"])),
    tolerance = 1e-12
  )
  expect_error(
    predict(fit.sc, newdata = data.frame(x = nd$x)),
    "must include columns"
  )

  bw.pl <- npplregbw(
    xdat = data.frame(z = z),
    ydat = y,
    zdat = data.frame(x = x),
    bws = matrix(c(0.25, 0.25), nrow = 2),
    bandwidth.compute = FALSE
  )
  fit.pl <- npplreg(
    bws = bw.pl,
    txdat = data.frame(z = z),
    tydat = y,
    tzdat = data.frame(x = x)
  )
  expect_equal(
    as.numeric(predict(fit.pl, newdata = nd)),
    as.numeric(predict(fit.pl, exdat = nd["z"], ezdat = nd["x"])),
    tolerance = 1e-12
  )
  expect_error(
    predict(fit.pl, newdata = data.frame(x = nd$x)),
    "must include columns"
  )
})
