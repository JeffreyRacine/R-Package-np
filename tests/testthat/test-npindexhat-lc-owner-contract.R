library(np)

test_that("npindexhat lc owner preserves fit and matrix/apply parity across bwtypes", {
  set.seed(20260315)
  n <- 90L
  x1 <- runif(n, -1, 1)
  x2 <- rnorm(n)
  y <- sin(2 * (0.7 * x1 - 0.3 * x2)) + 0.25 * x1 + rnorm(n, sd = 0.04)

  tx <- data.frame(x1 = x1, x2 = x2)
  ex <- data.frame(
    x1 = seq(min(x1) * 0.9, max(x1) * 0.9, length.out = 30L),
    x2 = seq(quantile(x2, 0.2), quantile(x2, 0.8), length.out = 30L)
  )

  for (bwtype in c("fixed", "generalized_nn", "adaptive_nn")) {
    bw <- npindexbw(
      xdat = tx,
      ydat = y,
      regtype = "lc",
      bwtype = bwtype,
      bandwidth.compute = FALSE,
      bws = c(1, 1, if (identical(bwtype, "fixed")) 0.22 else 9L)
    )

    fit.mean <- npindex(
      bws = bw,
      txdat = tx,
      tydat = y,
      exdat = ex,
      gradients = FALSE,
      errors = FALSE
    )
    fit.grad <- npindex(
      bws = bw,
      txdat = tx,
      tydat = y,
      exdat = ex,
      gradients = TRUE,
      errors = FALSE
    )

    apply.mean <- npindexhat(
      bws = bw,
      txdat = tx,
      exdat = ex,
      y = y,
      output = "apply",
      s = 0L
    )
    matrix.mean <- npindexhat(
      bws = bw,
      txdat = tx,
      exdat = ex,
      output = "matrix",
      s = 0L
    )
    apply.grad <- npindexhat(
      bws = bw,
      txdat = tx,
      exdat = ex,
      y = y,
      output = "apply",
      s = 1L
    )
    matrix.grad <- npindexhat(
      bws = bw,
      txdat = tx,
      exdat = ex,
      output = "matrix",
      s = 1L
    )

    expect_equal(as.vector(apply.mean), as.vector(fit.mean$mean), tolerance = 1e-8, info = paste("mean", bwtype))
    expect_equal(as.vector(matrix.mean %*% y), as.vector(fit.mean$mean), tolerance = 1e-8, info = paste("mean matrix", bwtype))
    expect_equal(as.vector(apply.mean), as.vector(matrix.mean %*% y), tolerance = 1e-10, info = paste("mean parity", bwtype))

    expect_equal(as.vector(apply.grad), as.vector(fit.grad$grad[, 1L]), tolerance = 1e-8, info = paste("grad", bwtype))
    expect_equal(as.vector(matrix.grad %*% y), as.vector(fit.grad$grad[, 1L]), tolerance = 1e-8, info = paste("grad matrix", bwtype))
    expect_equal(as.vector(apply.grad), as.vector(matrix.grad %*% y), tolerance = 1e-10, info = paste("grad parity", bwtype))
  }
})
