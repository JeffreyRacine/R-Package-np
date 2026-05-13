library(np)

test_that("semiparametric formula estimators reject malformed newdata before parent-frame lookup", {
  set.seed(20260513)
  dat <- data.frame(
    y = rnorm(60),
    x1 = runif(60),
    x2 = runif(60),
    z = runif(60)
  )

  bw.si <- npindexbw(y ~ x1 + x2, data = dat,
                     bws = c(0.25, 0.25, 1),
                     bandwidth.compute = FALSE)
  fit.si <- npindex(bws = bw.si, gradients = TRUE, errors = TRUE)
  expect_error(npindex(bws = bw.si, newdata = data.frame(q = 1:2)),
               "newdata must contain columns")
  expect_error(predict(fit.si, newdata = data.frame(q = 1:2)),
               "newdata must contain columns")

  bw.sc <- npscoefbw(y ~ x1 | z, data = dat, bws = 0.25,
                     bandwidth.compute = FALSE)
  fit.sc <- npscoef(bws = bw.sc, errors = TRUE, iterate = FALSE)
  expect_error(npscoef(bws = bw.sc, newdata = data.frame(q = 1:2),
                       errors = FALSE, iterate = FALSE),
               "newdata must contain columns")
  expect_error(predict(fit.sc, newdata = data.frame(q = 1:2)),
               "newdata must contain columns")

  bw.pl <- npplregbw(y ~ x2 | z, data = dat,
                     bws = matrix(c(0.25, 0.25), nrow = 2),
                     bandwidth.compute = FALSE)
  fit.pl <- npplreg(bws = bw.pl)
  expect_error(npplreg(bws = bw.pl, newdata = data.frame(z = c(0.3, 0.5))),
               "newdata must contain columns")
  expect_error(predict(fit.pl, newdata = data.frame(x2 = c(0.3, 0.5))),
               "newdata must contain columns")
})

test_that("semiparametric predict methods retain native evaluation argument precedence", {
  set.seed(20260513)
  dat <- data.frame(
    y = rnorm(60),
    x1 = runif(60),
    x2 = runif(60),
    z = runif(60)
  )
  nd.standard <- data.frame(x1 = c(0.2, 0.4, 0.6),
                            x2 = c(0.1, 0.3, 0.5),
                            z = c(0.7, 0.5, 0.3))
  nd.native <- data.frame(x1 = c(0.15, 0.35, 0.55),
                          x2 = c(0.12, 0.32, 0.52),
                          z = c(0.65, 0.45, 0.25))

  bw.si <- npindexbw(y ~ x1 + x2, data = dat,
                     bws = c(0.25, 0.25, 1),
                     bandwidth.compute = FALSE)
  fit.si <- npindex(bws = bw.si)
  expect_equal(
    as.numeric(predict(fit.si,
                       newdata = nd.standard[c("x1", "x2")],
                       exdat = nd.native[c("x1", "x2")])),
    as.numeric(predict(fit.si, exdat = nd.native[c("x1", "x2")])),
    tolerance = 1e-12
  )

  bw.sc <- npscoefbw(y ~ x1 | z, data = dat, bws = 0.25,
                     bandwidth.compute = FALSE)
  fit.sc <- npscoef(bws = bw.sc, errors = TRUE, iterate = FALSE)
  expect_equal(
    as.numeric(predict(fit.sc,
                       newdata = nd.standard[c("x1", "z")],
                       exdat = nd.native["x1"],
                       ezdat = nd.native["z"])),
    as.numeric(predict(fit.sc, exdat = nd.native["x1"], ezdat = nd.native["z"])),
    tolerance = 1e-12
  )

  bw.pl <- npplregbw(y ~ x2 | z, data = dat,
                     bws = matrix(c(0.25, 0.25), nrow = 2),
                     bandwidth.compute = FALSE)
  fit.pl <- npplreg(bws = bw.pl)
  expect_equal(
    as.numeric(predict(fit.pl,
                       newdata = nd.standard[c("x2", "z")],
                       exdat = nd.native["x2"],
                       ezdat = nd.native["z"])),
    as.numeric(predict(fit.pl, exdat = nd.native["x2"], ezdat = nd.native["z"])),
    tolerance = 1e-12
  )
})

test_that("semiparametric S3 logical arguments fail clearly", {
  set.seed(20260513)
  dat <- data.frame(
    y = rnorm(60),
    x1 = runif(60),
    x2 = runif(60),
    z = runif(60)
  )

  bw.si <- npindexbw(y ~ x1 + x2, data = dat,
                     bws = c(0.25, 0.25, 1),
                     bandwidth.compute = FALSE)
  fit.si <- npindex(bws = bw.si, gradients = TRUE, errors = TRUE)
  expect_error(predict(fit.si, se.fit = "yes"),
               "'se.fit' must be TRUE or FALSE", fixed = TRUE)
  expect_error(gradients(fit.si, errors = "yes"),
               "'errors' must be TRUE or FALSE", fixed = TRUE)

  bw.sc <- npscoefbw(y ~ x1 | z, data = dat, bws = 0.25,
                     bandwidth.compute = FALSE)
  fit.sc <- npscoef(bws = bw.sc, errors = TRUE, iterate = FALSE)
  expect_error(predict(fit.sc, se.fit = "yes"),
               "'se.fit' must be TRUE or FALSE", fixed = TRUE)

  bw.pl <- npplregbw(y ~ x2 | z, data = dat,
                     bws = matrix(c(0.25, 0.25), nrow = 2),
                     bandwidth.compute = FALSE)
  fit.pl <- npplreg(bws = bw.pl)
  expect_error(predict(fit.pl, se.fit = "yes"),
               "'se.fit' must be TRUE or FALSE", fixed = TRUE)
})
