test_that("npreghat reproduces npreg fitted values for mixed-data local constant", {
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

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
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

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
