test_that("adaptive ll matrix hat matches apply and npreg for mean and gradient", {
  set.seed(42)

  n <- 80
  x <- runif(n)
  y <- x + rnorm(n)
  tx <- data.frame(x = x)
  ex <- data.frame(x = seq(min(x), max(x), length.out = 25))

  bw <- npregbw(y ~ x, regtype = "ll", bwtype = "adaptive_nn")

  mean.apply <- npreghat(
    bws = bw,
    txdat = tx,
    exdat = ex,
    y = y,
    output = "apply"
  )
  mean.matrix <- npreghat(
    bws = bw,
    txdat = tx,
    exdat = ex,
    output = "matrix"
  )

  grad.apply <- npreghat(
    bws = bw,
    txdat = tx,
    exdat = ex,
    y = y,
    output = "apply",
    s = 1
  )
  grad.matrix <- npreghat(
    bws = bw,
    txdat = tx,
    exdat = ex,
    output = "matrix",
    s = 1
  )

  np.mean <- npreg(txdat = tx, tydat = y, bws = bw, exdat = ex)$mean
  np.grad <- npreg(txdat = tx, tydat = y, bws = bw, exdat = ex, gradients = TRUE)$grad[, 1L]

  expect_equal(drop(mean.matrix %*% y), mean.apply, tolerance = 1e-12)
  expect_equal(drop(mean.matrix %*% y), np.mean, tolerance = 1e-12)
  expect_equal(drop(grad.matrix %*% y), grad.apply, tolerance = 1e-12)
  expect_equal(drop(grad.matrix %*% y), np.grad, tolerance = 1e-12)
})
