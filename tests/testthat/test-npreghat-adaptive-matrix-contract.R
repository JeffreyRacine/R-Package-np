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

test_that("adaptive lc matrix hat matches fitted values and predict on train and eval data", {
  set.seed(20260315)

  n <- 64L
  x <- sort(runif(n))
  y <- sin(2 * pi * x) + 0.2 * x + rnorm(n, sd = 0.03)
  tx <- data.frame(x = x)
  ex <- data.frame(x = seq(0.05, 0.95, length.out = 21L))

  bw <- npregbw(y ~ x, regtype = "lc", bwtype = "adaptive_nn", nmulti = 1L)

  fit.train <- npreg(bws = bw, txdat = tx, tydat = y)
  fit.eval <- npreg(bws = bw, txdat = tx, tydat = y, exdat = ex)

  hat.train.matrix <- npreghat(bws = bw, txdat = tx, output = "matrix")
  hat.train.apply <- npreghat(bws = bw, txdat = tx, y = y, output = "apply")
  hat.eval.matrix <- npreghat(bws = bw, txdat = tx, exdat = ex, output = "matrix")
  hat.eval.apply <- npreghat(bws = bw, txdat = tx, exdat = ex, y = y, output = "apply")

  expect_equal(drop(hat.train.matrix %*% y), as.vector(fit.train$mean), tolerance = 1e-12)
  expect_equal(as.vector(hat.train.apply), as.vector(fit.train$mean), tolerance = 1e-12)
  expect_equal(drop(hat.train.matrix %*% y), as.vector(hat.train.apply), tolerance = 1e-12)

  expect_equal(drop(hat.eval.matrix %*% y), as.vector(fit.eval$mean), tolerance = 1e-12)
  expect_equal(as.vector(hat.eval.apply), as.vector(fit.eval$mean), tolerance = 1e-12)
  expect_equal(drop(hat.eval.matrix %*% y), as.vector(hat.eval.apply), tolerance = 1e-12)
  expect_equal(drop(hat.eval.matrix %*% y), as.vector(predict(fit.train, newdata = ex)), tolerance = 1e-12)
})

test_that("adaptive lp matrix hat matches apply and npreg on train and eval data across degrees", {
  set.seed(20260315)

  n <- 72L
  x <- sort(rnorm(n))
  y <- x^3 - 0.4 * x + rnorm(n, sd = 0.04)
  tx <- data.frame(x = x)
  ex <- data.frame(x = seq(min(x), max(x), length.out = 27L))

  for (degree in 1:3) {
    bw <- npregbw(
      xdat = tx,
      ydat = y,
      regtype = "lp",
      degree = degree,
      bwtype = "adaptive_nn",
      bandwidth.compute = FALSE,
      bws = 18L
    )

    fit.train <- npreg(bws = bw, txdat = tx, tydat = y, gradients = TRUE)
    fit.eval <- npreg(bws = bw, txdat = tx, tydat = y, exdat = ex, gradients = TRUE)

    mean.train.apply <- npreghat(bws = bw, txdat = tx, y = y, output = "apply")
    mean.train.matrix <- npreghat(bws = bw, txdat = tx, output = "matrix")
    mean.eval.apply <- npreghat(bws = bw, txdat = tx, exdat = ex, y = y, output = "apply")
    mean.eval.matrix <- npreghat(bws = bw, txdat = tx, exdat = ex, output = "matrix")

    grad.train.apply <- npreghat(bws = bw, txdat = tx, y = y, output = "apply", s = 1L)
    grad.train.matrix <- npreghat(bws = bw, txdat = tx, output = "matrix", s = 1L)
    grad.eval.apply <- npreghat(bws = bw, txdat = tx, exdat = ex, y = y, output = "apply", s = 1L)
    grad.eval.matrix <- npreghat(bws = bw, txdat = tx, exdat = ex, output = "matrix", s = 1L)

    expect_equal(drop(mean.train.matrix %*% y), as.vector(fit.train$mean), tolerance = 1e-10)
    expect_equal(as.vector(mean.train.apply), as.vector(fit.train$mean), tolerance = 1e-10)
    expect_equal(drop(mean.train.matrix %*% y), as.vector(mean.train.apply), tolerance = 1e-10)

    expect_equal(drop(mean.eval.matrix %*% y), as.vector(fit.eval$mean), tolerance = 1e-10)
    expect_equal(as.vector(mean.eval.apply), as.vector(fit.eval$mean), tolerance = 1e-10)
    expect_equal(drop(mean.eval.matrix %*% y), as.vector(mean.eval.apply), tolerance = 1e-10)

    expect_equal(drop(grad.train.matrix %*% y), as.vector(fit.train$grad[, 1L]), tolerance = 1e-8)
    expect_equal(as.vector(grad.train.apply), as.vector(fit.train$grad[, 1L]), tolerance = 1e-8)
    expect_equal(drop(grad.train.matrix %*% y), as.vector(grad.train.apply), tolerance = 1e-8)

    expect_equal(drop(grad.eval.matrix %*% y), as.vector(fit.eval$grad[, 1L]), tolerance = 1e-8)
    expect_equal(as.vector(grad.eval.apply), as.vector(fit.eval$grad[, 1L]), tolerance = 1e-8)
    expect_equal(drop(grad.eval.matrix %*% y), as.vector(grad.eval.apply), tolerance = 1e-8)
  }
})

test_that("fixed and generalized lc/ll matrix hats match apply and npreg on train and eval data", {
  set.seed(20260316)

  n <- 96L
  x <- sort(rnorm(n))
  y <- x^2 + 0.15 * x + rnorm(n, sd = 0.05)
  tx <- data.frame(x = x)
  ex <- data.frame(x = seq(min(x), max(x), length.out = 31L))

  for (regtype in c("lc", "ll")) {
    for (bwtype in c("fixed", "generalized_nn")) {
      bw <- npregbw(
        xdat = tx,
        ydat = y,
        regtype = regtype,
        bwtype = bwtype,
        bandwidth.compute = FALSE,
        bws = if (identical(bwtype, "fixed")) 0.55 else 9L
      )

      fit.train <- npreg(bws = bw, txdat = tx, tydat = y, gradients = identical(regtype, "ll"))
      fit.eval <- npreg(bws = bw, txdat = tx, tydat = y, exdat = ex, gradients = identical(regtype, "ll"))

      mean.train.matrix <- npreghat(bws = bw, txdat = tx, output = "matrix")
      mean.train.apply <- npreghat(bws = bw, txdat = tx, y = y, output = "apply")
      mean.eval.matrix <- npreghat(bws = bw, txdat = tx, exdat = ex, output = "matrix")
      mean.eval.apply <- npreghat(bws = bw, txdat = tx, exdat = ex, y = y, output = "apply")

      expect_equal(drop(mean.train.matrix %*% y), as.vector(fit.train$mean), tolerance = 1e-10)
      expect_equal(as.vector(mean.train.apply), as.vector(fit.train$mean), tolerance = 1e-10)
      expect_equal(drop(mean.eval.matrix %*% y), as.vector(fit.eval$mean), tolerance = 1e-10)
      expect_equal(as.vector(mean.eval.apply), as.vector(fit.eval$mean), tolerance = 1e-10)

      if (identical(regtype, "ll")) {
        grad.train.matrix <- npreghat(bws = bw, txdat = tx, output = "matrix", s = 1L)
        grad.train.apply <- npreghat(bws = bw, txdat = tx, y = y, output = "apply", s = 1L)
        grad.eval.matrix <- npreghat(bws = bw, txdat = tx, exdat = ex, output = "matrix", s = 1L)
        grad.eval.apply <- npreghat(bws = bw, txdat = tx, exdat = ex, y = y, output = "apply", s = 1L)

        expect_equal(drop(grad.train.matrix %*% y), as.vector(fit.train$grad[, 1L]), tolerance = 1e-8)
        expect_equal(as.vector(grad.train.apply), as.vector(fit.train$grad[, 1L]), tolerance = 1e-8)
        expect_equal(drop(grad.eval.matrix %*% y), as.vector(fit.eval$grad[, 1L]), tolerance = 1e-8)
        expect_equal(as.vector(grad.eval.apply), as.vector(fit.eval$grad[, 1L]), tolerance = 1e-8)
      }
    }
  }
})

test_that("fixed lp matrix hat matches apply and npreg across higher degrees", {
  set.seed(20260317)

  n <- 96L
  x <- sort(rnorm(n))
  y <- x^3 - 0.4 * x + rnorm(n, sd = 0.05)
  tx <- data.frame(x = x)
  ex <- data.frame(x = seq(min(x), max(x), length.out = 31L))

  for (degree in 2:3) {
    bw <- npregbw(
      xdat = tx,
      ydat = y,
      regtype = "lp",
      degree = degree,
      basis = "glp",
      bernstein.basis = FALSE,
      bwtype = "fixed",
      bandwidth.compute = FALSE,
      bws = 0.55
    )

    fit.train <- npreg(bws = bw, txdat = tx, tydat = y, gradients = TRUE)
    fit.eval <- npreg(bws = bw, txdat = tx, tydat = y, exdat = ex, gradients = TRUE)

    mean.train.matrix <- npreghat(bws = bw, txdat = tx, output = "matrix")
    mean.train.apply <- npreghat(bws = bw, txdat = tx, y = y, output = "apply")
    mean.eval.matrix <- npreghat(bws = bw, txdat = tx, exdat = ex, output = "matrix")
    mean.eval.apply <- npreghat(bws = bw, txdat = tx, exdat = ex, y = y, output = "apply")

    grad.train.matrix <- npreghat(bws = bw, txdat = tx, output = "matrix", s = 1L)
    grad.train.apply <- npreghat(bws = bw, txdat = tx, y = y, output = "apply", s = 1L)
    grad.eval.matrix <- npreghat(bws = bw, txdat = tx, exdat = ex, output = "matrix", s = 1L)
    grad.eval.apply <- npreghat(bws = bw, txdat = tx, exdat = ex, y = y, output = "apply", s = 1L)

    expect_equal(drop(mean.train.matrix %*% y), as.vector(fit.train$mean), tolerance = 1e-10)
    expect_equal(as.vector(mean.train.apply), as.vector(fit.train$mean), tolerance = 1e-10)
    expect_equal(drop(mean.eval.matrix %*% y), as.vector(fit.eval$mean), tolerance = 1e-10)
    expect_equal(as.vector(mean.eval.apply), as.vector(fit.eval$mean), tolerance = 1e-10)

    expect_equal(drop(grad.train.matrix %*% y), as.vector(fit.train$grad[, 1L]), tolerance = 1e-8)
    expect_equal(as.vector(grad.train.apply), as.vector(fit.train$grad[, 1L]), tolerance = 1e-8)
    expect_equal(drop(grad.eval.matrix %*% y), as.vector(fit.eval$grad[, 1L]), tolerance = 1e-8)
    expect_equal(as.vector(grad.eval.apply), as.vector(fit.eval$grad[, 1L]), tolerance = 1e-8)
  }
})
