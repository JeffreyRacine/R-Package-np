library(np)

test_that("core formula estimators reject malformed newdata before parent-frame lookup", {
  set.seed(20260513)
  dat <- data.frame(x = runif(45), y = runif(45))
  bad <- data.frame(z = dat$x[1:3])

  bw.reg <- npregbw(y ~ x, data = dat, bws = 0.25, bandwidth.compute = FALSE)
  fit.reg <- npreg(bws = bw.reg, gradients = TRUE)
  expect_error(npreg(bws = bw.reg, newdata = bad), "newdata must contain columns")
  expect_error(predict(fit.reg, newdata = bad), "newdata must contain columns")

  bw.den <- npudensbw(~ x, data = dat, bws = 0.25, bandwidth.compute = FALSE)
  fit.den <- npudens(bws = bw.den)
  expect_error(npudens(bws = bw.den, newdata = bad), "newdata must contain columns")
  expect_error(predict(fit.den, newdata = bad), "newdata must contain columns")

  bw.dist <- npudistbw(~ x, data = dat, bws = 0.25, bandwidth.compute = FALSE)
  fit.dist <- npudist(bws = bw.dist)
  expect_error(npudist(bws = bw.dist, newdata = bad), "newdata must contain columns")
  expect_error(predict(fit.dist, newdata = bad), "newdata must contain columns")

  bw.cd <- npcdensbw(y ~ x, data = dat, bws = c(0.25, 0.25), bandwidth.compute = FALSE)
  fit.cd <- npcdens(bws = bw.cd, gradients = TRUE)
  expect_error(npcdens(bws = bw.cd, newdata = bad), "newdata must contain columns")
  expect_error(predict(fit.cd, newdata = bad), "newdata must contain columns")

  bw.cdist <- npcdistbw(y ~ x, data = dat, bws = c(0.25, 0.25), bandwidth.compute = FALSE)
  fit.cdist <- npcdist(bws = bw.cdist, gradients = TRUE)
  expect_error(npcdist(bws = bw.cdist, newdata = bad), "newdata must contain columns")
  expect_error(predict(fit.cdist, newdata = bad), "newdata must contain columns")
})

test_that("core predict methods retain native evaluation argument precedence", {
  set.seed(20260513)
  dat <- data.frame(x = runif(45), y = runif(45))
  nd.standard <- data.frame(x = c(0.2, 0.5, 0.8))
  nd.native <- data.frame(x = c(0.15, 0.35, 0.65))
  cd.standard <- data.frame(y = c(0.2, 0.4, 0.7), x = c(0.2, 0.5, 0.8))
  cd.native <- data.frame(y = c(0.1, 0.3, 0.6), x = c(0.15, 0.35, 0.65))

  bw.reg <- npregbw(y ~ x, data = dat, bws = 0.25, bandwidth.compute = FALSE)
  fit.reg <- npreg(bws = bw.reg)
  expect_equal(
    as.numeric(predict(fit.reg, newdata = nd.standard, exdat = nd.native)),
    as.numeric(predict(fit.reg, exdat = nd.native)),
    tolerance = 1e-12
  )

  bw.den <- npudensbw(~ x, data = dat, bws = 0.25, bandwidth.compute = FALSE)
  fit.den <- npudens(bws = bw.den)
  expect_equal(
    as.numeric(predict(fit.den, newdata = nd.standard, edat = nd.native)),
    as.numeric(predict(fit.den, edat = nd.native)),
    tolerance = 1e-12
  )

  bw.dist <- npudistbw(~ x, data = dat, bws = 0.25, bandwidth.compute = FALSE)
  fit.dist <- npudist(bws = bw.dist)
  expect_equal(
    as.numeric(predict(fit.dist, newdata = nd.standard, edat = nd.native)),
    as.numeric(predict(fit.dist, edat = nd.native)),
    tolerance = 1e-12
  )

  bw.cd <- npcdensbw(y ~ x, data = dat, bws = c(0.25, 0.25), bandwidth.compute = FALSE)
  fit.cd <- npcdens(bws = bw.cd)
  expect_equal(
    as.numeric(predict(fit.cd, newdata = cd.standard, eydat = cd.native["y"], exdat = cd.native["x"])),
    as.numeric(predict(fit.cd, eydat = cd.native["y"], exdat = cd.native["x"])),
    tolerance = 1e-12
  )

  bw.cdist <- npcdistbw(y ~ x, data = dat, bws = c(0.25, 0.25), bandwidth.compute = FALSE)
  fit.cdist <- npcdist(bws = bw.cdist)
  expect_equal(
    as.numeric(predict(fit.cdist, newdata = cd.standard, eydat = cd.native["y"], exdat = cd.native["x"])),
    as.numeric(predict(fit.cdist, eydat = cd.native["y"], exdat = cd.native["x"])),
    tolerance = 1e-12
  )
})

test_that("core S3 logical arguments fail clearly", {
  set.seed(20260513)
  dat <- data.frame(x = runif(45), y = runif(45))

  bw.reg <- npregbw(y ~ x, data = dat, bws = 0.25, bandwidth.compute = FALSE)
  fit.reg <- npreg(bws = bw.reg, gradients = TRUE)
  expect_error(predict(fit.reg, se.fit = "yes"), "'se.fit' must be TRUE or FALSE", fixed = TRUE)
  expect_error(gradients(fit.reg, errors = "yes"), "'errors' must be TRUE or FALSE", fixed = TRUE)

  bw.den <- npudensbw(~ x, data = dat, bws = 0.25, bandwidth.compute = FALSE)
  fit.den <- npudens(bws = bw.den)
  expect_error(predict(fit.den, se.fit = "yes"), "'se.fit' must be TRUE or FALSE", fixed = TRUE)

  bw.dist <- npudistbw(~ x, data = dat, bws = 0.25, bandwidth.compute = FALSE)
  fit.dist <- npudist(bws = bw.dist)
  expect_error(predict(fit.dist, se.fit = "yes"), "'se.fit' must be TRUE or FALSE", fixed = TRUE)

  bw.cd <- npcdensbw(y ~ x, data = dat, bws = c(0.25, 0.25), bandwidth.compute = FALSE)
  fit.cd <- npcdens(bws = bw.cd, gradients = TRUE)
  expect_error(predict(fit.cd, se.fit = "yes"), "'se.fit' must be TRUE or FALSE", fixed = TRUE)
  expect_error(gradients(fit.cd, errors = "yes"), "'errors' must be TRUE or FALSE", fixed = TRUE)

  bw.cdist <- npcdistbw(y ~ x, data = dat, bws = c(0.25, 0.25), bandwidth.compute = FALSE)
  fit.cdist <- npcdist(bws = bw.cdist, gradients = TRUE)
  expect_error(predict(fit.cdist, se.fit = "yes"), "'se.fit' must be TRUE or FALSE", fixed = TRUE)
  expect_error(gradients(fit.cdist, errors = "yes"), "'errors' must be TRUE or FALSE", fixed = TRUE)
})
