test_that("npreghat reproduces npreg fitted values for mixed-data local constant", {
  set.seed(20260223)
  n <- 120
  x <- runif(n)
  u <- factor(sample(c("a", "b", "c"), n, replace = TRUE))
  o <- ordered(sample(1:3, n, replace = TRUE))
  y <- sin(2 * pi * x) + as.numeric(u) / 3 + as.numeric(o) / 5 + rnorm(n, sd = 0.05)

  tx <- data.frame(x = x, u = u, o = o)
  bw <- npregbw(
    xdat = tx,
    ydat = y,
    bws = c(0.2, 0.4, 0.4),
    regtype = "lc",
    bandwidth.compute = FALSE
  )

  fit <- npreg(txdat = tx, tydat = y, bws = bw, warn.glp.gradient = FALSE)
  H <- npreghat(bws = bw, txdat = tx)

  expect_s3_class(H, "npreghat")
  expect_equal(as.vector(H %*% y), as.vector(fit$mean), tolerance = 1e-8)
})

test_that("npreghat supports lp/ll derivatives and matrix apply mode", {
  set.seed(777)
  n <- 150
  x <- sort(runif(n))
  y <- cos(2 * pi * x) + rnorm(n, sd = 0.03)
  tx <- data.frame(x = x)
  ex <- data.frame(x = seq(min(x), max(x), length.out = 40))

  bw <- npregbw(
    xdat = tx,
    ydat = y,
    bws = 0.2,
    regtype = "ll",
    bandwidth.compute = FALSE
  )

  fit <- npreg(
    txdat = tx,
    tydat = y,
    exdat = ex,
    bws = bw,
    gradients = TRUE,
    gradient.order = 1L,
    warn.glp.gradient = FALSE
  )

  H0 <- npreghat(bws = bw, txdat = tx, exdat = ex)
  H1 <- npreghat(bws = bw, txdat = tx, exdat = ex, s = 1L)

  expect_equal(as.vector(H0 %*% y), as.vector(fit$mean), tolerance = 1e-8)
  expect_equal(as.vector(H1 %*% y), as.vector(fit$grad[, 1]), tolerance = 1e-6)

  yboot <- cbind(y, y + 0.1)
  hy.apply <- npreghat(
    bws = bw,
    txdat = tx,
    exdat = ex,
    y = yboot,
    output = "apply"
  )

  expect_true(isTRUE(all.equal(
    hy.apply,
    H0 %*% yboot,
    tolerance = 1e-10,
    check.attributes = FALSE
  )))
})

test_that("npreghat leave.one.out is honored and predict reuses it safely", {
  set.seed(90210)
  n <- 80
  x <- runif(n)
  y <- sin(2 * pi * x) + rnorm(n, sd = 0.02)
  tx <- data.frame(x = x)
  bw <- npregbw(
    xdat = tx,
    ydat = y,
    bws = 0.18,
    regtype = "lc",
    bandwidth.compute = FALSE
  )

  H.in <- npreghat(bws = bw, txdat = tx, leave.one.out = FALSE)
  H.loo <- npreghat(bws = bw, txdat = tx, leave.one.out = TRUE)

  expect_gt(max(abs(H.in - H.loo)), 1e-8)
  expect_lt(max(abs(diag(H.loo))), 1e-12)
  expect_error(
    npreghat(bws = bw, txdat = tx, exdat = tx[1:10, , drop = FALSE], leave.one.out = TRUE),
    "you may not specify 'leave.one.out = TRUE' and provide evaluation data"
  )

  hy <- predict(H.loo, y = y, output = "apply")
  expect_equal(as.vector(hy), as.vector(H.loo %*% y), tolerance = 1e-10)
})

test_that("npreghat lp bernstein path matches predict semantics", {
  set.seed(20260305)
  n <- 220
  x <- sort(runif(n))
  y <- sin(2 * pi * x) + rnorm(n, sd = 0.08)
  tx <- data.frame(x = x)
  ex <- data.frame(x = seq(min(x), max(x), length.out = 60))

  bw.raw <- npregbw(
    xdat = tx,
    ydat = y,
    bws = 0.18,
    regtype = "lp",
    degree = 3L,
    bernstein = FALSE,
    bandwidth.compute = FALSE
  )
  bw.bern <- npregbw(
    xdat = tx,
    ydat = y,
    bws = 0.18,
    regtype = "lp",
    degree = 3L,
    bernstein = TRUE,
    bandwidth.compute = FALSE
  )

  fit.raw <- npreg(txdat = tx, tydat = y, exdat = ex, bws = bw.raw, gradients = FALSE, warn.glp.gradient = FALSE)
  fit.bern <- npreg(txdat = tx, tydat = y, exdat = ex, bws = bw.bern, gradients = FALSE, warn.glp.gradient = FALSE)

  hat.raw <- npreghat(bws = bw.raw, txdat = tx, exdat = ex, y = y, output = "apply", s = 0L)
  hat.bern <- npreghat(bws = bw.bern, txdat = tx, exdat = ex, y = y, output = "apply", s = 0L)

  expect_equal(as.vector(hat.raw), as.vector(fit.raw$mean), tolerance = 1e-8)
  expect_equal(as.vector(hat.bern), as.vector(fit.bern$mean), tolerance = 1e-8)
  expect_equal(as.vector(fit.bern$mean), as.vector(fit.raw$mean), tolerance = 1e-8)
})
