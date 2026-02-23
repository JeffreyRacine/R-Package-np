test_that("npreg returns higher-order lp gradients via hat-operator fallback", {
  set.seed(99)
  n <- 140
  x <- sort(runif(n))
  y <- sin(2 * pi * x) + rnorm(n, sd = 0.02)
  tx <- data.frame(x = x)
  ex <- data.frame(x = seq(min(x), max(x), length.out = 35))

  bw <- npregbw(
    xdat = tx,
    ydat = y,
    regtype = "lp",
    degree = 2L,
    basis = "glp",
    bws = 0.25,
    bandwidth.compute = FALSE
  )

  fit <- suppressWarnings(npreg(
    txdat = tx,
    tydat = y,
    exdat = ex,
    bws = bw,
    gradients = TRUE,
    gradient.order = 2L
  ))

  H2 <- npreghat(bws = bw, txdat = tx, exdat = ex, s = 2L)
  g2 <- as.vector(H2 %*% y)

  expect_false(all(is.na(fit$grad[, 1])))
  expect_equal(as.vector(fit$grad[, 1]), g2, tolerance = 1e-6)
})
