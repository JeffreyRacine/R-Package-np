test_that("Klein-Spady covariance uses the centered information score", {
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  set.seed(6106)
  n <- 160L
  x <- data.frame(x = rnorm(n), z = rbinom(n, 1L, .5), w = rnorm(n))
  y <- rbinom(n, 1L, plogis(x$x + .6*x$z - .2*x$w))
  h <- .65
  b <- npindexbw(xdat = x, ydat = y, bws = c(1, .6, -.2, h),
                 method = "kleinspady", bandwidth.compute = FALSE)
  f <- npindex(bws = b, txdat = x, tydat = y, gradients = TRUE)
  v <- drop(as.matrix(x) %*% b$beta)
  # Independent fixed-Gaussian conditional means; no package moment helper.
  K <- dnorm(outer(v, v, `-`) / h)
  Z <- as.matrix(x[-1L])
  mu <- sweep(crossprod(K, Z), 1L, colSums(K), `/`)
  score <- (Z - mu) * gradients(f)[, 1L]
  p <- fitted(f)
  reference <- solve(crossprod(score, score / (p * (1 - p))))
  expect_equal(unname(vcov(f)[-1L, -1L]), unname(reference), tolerance = 1e-10)
  expect_identical(unname(vcov(f)[1L, ]), rep(0, 3L))
  shifted <- transform(x, z = z + 1, w = w + 2)
  g <- npindex(bws = b, txdat = shifted, tydat = y, gradients = TRUE)
  expect_equal(fitted(g), fitted(f), tolerance = 1e-12)
  expect_equal(vcov(g), vcov(f), tolerance = 1e-10)
  plain <- npindex(bws = b, txdat = x, tydat = y, se = FALSE, gradients = TRUE)
  expect_identical(fitted(plain), fitted(f))
  expect_null(plain$betavcov)
})
