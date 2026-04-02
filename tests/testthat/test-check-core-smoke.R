test_that("npudens core smoke stays alive", {
  set.seed(1)
  x <- data.frame(x = seq(0.1, 0.9, length.out = 24L))

  bw <- npudensbw(dat = x, bws = 0.2, bandwidth.compute = FALSE)
  fit <- npudens(bws = bw)

  expect_s3_class(fit, "npdensity")
  expect_equal(length(predict(fit)), nrow(x))
  expect_true(all(is.finite(predict(fit))))
})

test_that("npudist core smoke stays alive", {
  set.seed(2)
  x <- data.frame(x = seq(0.1, 0.9, length.out = 24L))

  bw <- npudistbw(dat = x, bws = 0.2, bandwidth.compute = FALSE)
  fit <- npudist(bws = bw)

  expect_s3_class(fit, "npdistribution")
  expect_equal(length(predict(fit)), nrow(x))
  expect_true(all(predict(fit) >= 0 & predict(fit) <= 1))
})

test_that("npreg core smoke stays alive", {
  set.seed(3)
  x <- data.frame(x = seq(0.1, 1.0, length.out = 28L))
  y <- sin(2 * pi * x$x)

  bw <- npregbw(xdat = x, ydat = y, bws = 0.25, bandwidth.compute = FALSE)
  fit <- npreg(bws = bw)

  expect_s3_class(fit, "npregression")
  expect_equal(length(predict(fit)), nrow(x))
  expect_true(all(is.finite(predict(fit))))
})

test_that("npcdens core smoke stays alive", {
  set.seed(4)
  x <- data.frame(x = seq(0.1, 1.0, length.out = 24L))
  y <- data.frame(y = x$x^2)

  bw <- npcdensbw(xdat = x, ydat = y, bws = c(0.25, 0.25), bandwidth.compute = FALSE)
  fit <- npcdens(bws = bw)

  expect_s3_class(fit, "condensity")
  expect_equal(length(predict(fit)), nrow(x))
  expect_true(all(is.finite(predict(fit))))
})

test_that("npcdist core smoke stays alive", {
  set.seed(5)
  x <- data.frame(x = seq(0.1, 1.0, length.out = 24L))
  y <- data.frame(y = x$x^2)

  bw <- npcdistbw(xdat = x, ydat = y, bws = c(0.25, 0.25), bandwidth.compute = FALSE)
  fit <- npcdist(bws = bw)

  expect_s3_class(fit, "condistribution")
  expect_equal(length(predict(fit)), nrow(x))
  expect_true(all(predict(fit) >= 0 & predict(fit) <= 1))
})

test_that("npplreg core smoke stays alive", {
  set.seed(6)
  n <- 24L
  z <- data.frame(z = seq(0.1, 1.0, length.out = n))
  x <- data.frame(x = seq(0.2, 1.1, length.out = n))
  y <- z$z^2 + 2 * x$x

  bw <- npplregbw(
    xdat = x,
    zdat = z,
    ydat = y,
    bws = matrix(c(0.25, 0.25), nrow = 2L),
    bandwidth.compute = FALSE
  )
  fit <- npplreg(bws = bw)

  expect_s3_class(fit, "plregression")
  expect_equal(length(predict(fit)), n)
  expect_true(all(is.finite(predict(fit))))
})

test_that("npindex core smoke stays alive", {
  set.seed(7)
  n <- 24L
  x <- data.frame(
    x1 = seq(0.1, 1.0, length.out = n),
    x2 = seq(1.0, 0.1, length.out = n)
  )
  y <- x$x1 - x$x2

  bw <- npindexbw(
    xdat = x,
    ydat = y,
    method = "ichimura",
    bws = c(1, 0.25, 0.25),
    bandwidth.compute = FALSE
  )
  fit <- npindex(bws = bw)

  expect_s3_class(fit, "singleindex")
  expect_equal(length(predict(fit)), n)
  expect_true(all(is.finite(predict(fit))))
})

test_that("npscoef core smoke stays alive", {
  set.seed(8)
  n <- 24L
  x <- data.frame(x = seq(0.1, 1.0, length.out = n))
  z <- data.frame(z = seq(1.0, 0.1, length.out = n))
  y <- x$x * z$z

  bw <- npscoefbw(
    xdat = x,
    zdat = z,
    ydat = y,
    bws = 0.25,
    bandwidth.compute = FALSE
  )
  fit <- npscoef(bws = bw, iterate = FALSE, errors = FALSE)

  expect_s3_class(fit, "smoothcoefficient")
  expect_equal(length(predict(fit)), n)
  expect_true(all(is.finite(predict(fit))))
})
