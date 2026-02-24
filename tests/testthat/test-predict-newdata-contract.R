library(np)

test_that("predict aliases newdata to native eval args for default npreg/npudens/npudist", {
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
