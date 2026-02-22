test_that("npksum numeric and formula interfaces agree", {
  set.seed(1)
  n <- 30
  x <- runif(n)
  y <- rnorm(n)
  dat <- data.frame(y = y, x = x)

  k1 <- npksum(txdat = x, tydat = y, bws = 0.5)
  k2 <- npksum(y ~ x, data = dat, bws = 0.5)

  expect_s3_class(k1, "npkernelsum")
  expect_s3_class(k2, "npkernelsum")
  expect_identical(k1$ksum, k2$ksum)
})

test_that("npksum preserves 1D exdat column naming", {
  set.seed(1)
  n <- 10
  x <- runif(n)
  ex <- seq(0, 1, length.out = 3)

  k <- npksum(txdat = x, exdat = ex, bws = 0.5)
  expect_identical(colnames(k$eval), "exdat")
})

test_that("npksum positional bws dispatch matches named bws", {
  set.seed(42)
  n <- 20
  x <- runif(n)
  y <- rnorm(n)

  k_named <- npksum(txdat = x, tydat = y, bws = 0.4)
  k_pos <- npksum(0.4, txdat = x, tydat = y)

  expect_s3_class(k_pos, "npkernelsum")
  expect_identical(k_named$ksum, k_pos$ksum)
})

test_that("npksum formula subset and na.action match explicit data path", {
  set.seed(7)
  n <- 24
  dat <- data.frame(x = runif(n), y = rnorm(n))
  dat$y[c(3, 12)] <- NA_real_

  k_formula <- npksum(y ~ x, data = dat, subset = x > 0.2, na.action = na.omit, bws = 0.45)

  dat2 <- na.omit(subset(dat, x > 0.2))
  k_explicit <- npksum(txdat = dat2$x, tydat = dat2$y, bws = 0.45)

  expect_identical(k_formula$ksum, k_explicit$ksum)
})
