test_that("formula npudensbw matches default interface with subset/na.action", {
  set.seed(20260222)
  dat <- data.frame(
    x1 = runif(30),
    x2 = rnorm(30)
  )
  dat$x2[c(4, 17)] <- NA_real_
  bw_formula <- np::npudensbw(
    ~ x1 + x2,
    data = dat,
    subset = x1 > 0.2,
    na.action = na.omit,
    bws = c(0.4, 0.6),
    bandwidth.compute = FALSE
  )

  mf <- model.frame(~ x1 + x2, data = dat, subset = x1 > 0.2, na.action = na.omit)
  bw_default <- np::npudensbw(
    dat = mf,
    bws = c(0.4, 0.6),
    bandwidth.compute = FALSE
  )

  expect_equal(as.numeric(bw_formula$bw), as.numeric(bw_default$bw))
})

test_that("formula npudistbw matches default interface with subset/na.action", {
  set.seed(20260222)
  dat <- data.frame(
    x1 = runif(30),
    x2 = rnorm(30)
  )
  dat$x1[c(2, 11)] <- NA_real_
  bw_formula <- np::npudistbw(
    ~ x1 + x2,
    data = dat,
    subset = x2 > -0.5,
    na.action = na.omit,
    bws = c(0.3, 0.5),
    bandwidth.compute = FALSE
  )

  mf <- model.frame(~ x1 + x2, data = dat, subset = x2 > -0.5, na.action = na.omit)
  bw_default <- np::npudistbw(
    dat = mf,
    bws = c(0.3, 0.5),
    bandwidth.compute = FALSE
  )

  expect_equal(as.numeric(bw_formula$bw), as.numeric(bw_default$bw))
})

test_that("formula npudistbw gdata path matches default interface", {
  set.seed(20260222)
  dat <- data.frame(
    x1 = runif(30),
    x2 = rnorm(30)
  )
  gdat <- data.frame(
    x1 = runif(22),
    x2 = rnorm(22)
  )

  bw_formula <- np::npudistbw(
    ~ x1 + x2,
    data = dat,
    gdata = gdat,
    bws = c(0.3, 0.5),
    bandwidth.compute = FALSE
  )

  mf <- model.frame(~ x1 + x2, data = dat)
  gmf <- model.frame(~ x1 + x2, data = gdat)
  bw_default <- np::npudistbw(
    dat = mf,
    gdat = gmf,
    bws = c(0.3, 0.5),
    bandwidth.compute = FALSE
  )

  expect_equal(as.numeric(bw_formula$bw), as.numeric(bw_default$bw))
})

test_that("formula npregbw matches default interface with subset/na.action", {
  set.seed(20260222)
  dat <- data.frame(
    y = rnorm(40),
    x1 = runif(40),
    x2 = rnorm(40)
  )
  dat$y[c(3, 9)] <- NA_real_
  bw_formula <- np::npregbw(
    y ~ x1 + x2,
    data = dat,
    subset = x1 < 0.85,
    na.action = na.omit,
    bws = c(0.45, 0.55),
    bandwidth.compute = FALSE,
    regtype = "ll"
  )

  mf <- model.frame(y ~ x1 + x2, data = dat, subset = x1 < 0.85, na.action = na.omit)
  bw_default <- np::npregbw(
    xdat = mf[, c("x1", "x2"), drop = FALSE],
    ydat = model.response(mf),
    bws = c(0.45, 0.55),
    bandwidth.compute = FALSE,
    regtype = "ll"
  )

  expect_equal(as.numeric(bw_formula$bw), as.numeric(bw_default$bw))
})

test_that("formula npcdistbw gdata path matches default interface", {
  set.seed(20260222)
  dat <- data.frame(
    y = rnorm(36),
    x = runif(36)
  )
  gdat <- data.frame(
    y = rnorm(18),
    x = runif(18)
  )

  bw_formula <- np::npcdistbw(
    y ~ x,
    data = dat,
    gdata = gdat,
    bws = c(0.35, 0.45),
    bandwidth.compute = FALSE
  )

  mf <- model.frame(y ~ x, data = dat)
  gmf <- model.frame(y ~ x, data = gdat)
  bw_default <- np::npcdistbw(
    xdat = mf[, "x", drop = FALSE],
    ydat = mf[, "y", drop = FALSE],
    gydat = gmf[, "y", drop = FALSE],
    bws = c(0.35, 0.45),
    bandwidth.compute = FALSE
  )

  expect_equal(as.numeric(bw_formula$bw), as.numeric(bw_default$bw))
})

test_that("formula npcdensbw matches default interface with subset/na.action", {
  set.seed(20260222)
  dat <- data.frame(
    y = rnorm(36),
    x = runif(36)
  )
  dat$y[c(5, 19)] <- NA_real_

  bw_formula <- np::npcdensbw(
    y ~ x,
    data = dat,
    subset = x > 0.25,
    na.action = na.omit,
    bws = c(0.35, 0.45),
    bandwidth.compute = FALSE
  )

  mf <- model.frame(y ~ x, data = dat, subset = x > 0.25, na.action = na.omit)
  bw_default <- np::npcdensbw(
    xdat = mf[, "x", drop = FALSE],
    ydat = mf[, "y", drop = FALSE],
    bws = c(0.35, 0.45),
    bandwidth.compute = FALSE
  )

  expect_equal(as.numeric(bw_formula$bw), as.numeric(bw_default$bw))
})

test_that("formula npindexbw matches default interface with subset/na.action", {
  set.seed(20260222)
  dat <- data.frame(
    y = rnorm(42),
    x1 = runif(42),
    x2 = rnorm(42)
  )
  dat$y[c(4, 21)] <- NA_real_

  bw_formula <- np::npindexbw(
    y ~ x1 + x2,
    data = dat,
    subset = x1 < 0.9,
    na.action = na.omit,
    bws = c(0.2, 0.3, 0.4),
    bandwidth.compute = FALSE,
    method = "ichimura",
    nmulti = 1
  )

  mf <- model.frame(y ~ x1 + x2, data = dat, subset = x1 < 0.9, na.action = na.omit)
  bw_default <- np::npindexbw(
    xdat = mf[, c("x1", "x2"), drop = FALSE],
    ydat = model.response(mf),
    bws = c(0.2, 0.3, 0.4),
    bandwidth.compute = FALSE,
    method = "ichimura",
    nmulti = 1
  )

  expect_equal(as.numeric(bw_formula$bw), as.numeric(bw_default$bw))
})
