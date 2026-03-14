library(np)

test_that("npindexhat ll owner preserves fit and matrix/apply parity across bwtypes", {
  set.seed(20260314)
  n <- 70
  x1 <- runif(n)
  x2 <- runif(n)
  y <- x1 + 0.5 * x2 + rnorm(n, sd = 0.04)

  tx <- data.frame(x1 = x1, x2 = x2)
  ex <- tx[seq_len(18), , drop = FALSE]

  for (bwtype in c("fixed", "generalized_nn", "adaptive_nn")) {
    bw <- npindexbw(
      xdat = tx,
      ydat = y,
      regtype = "ll",
      bwtype = bwtype,
      bandwidth.compute = FALSE,
      bws = c(1, 1, if (identical(bwtype, "fixed")) 0.18 else 9)
    )

    fit.mean <- npindex(
      bws = bw,
      txdat = tx,
      tydat = y,
      exdat = ex,
      gradients = FALSE
    )
    fit.grad <- npindex(
      bws = bw,
      txdat = tx,
      tydat = y,
      exdat = ex,
      gradients = TRUE
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

    expect_equal(as.vector(apply.grad), as.vector(fit.grad$grad[, 1]), tolerance = 1e-8, info = paste("grad", bwtype))
    expect_equal(as.vector(matrix.grad %*% y), as.vector(fit.grad$grad[, 1]), tolerance = 1e-8, info = paste("grad matrix", bwtype))
    expect_equal(as.vector(apply.grad), as.vector(matrix.grad %*% y), tolerance = 1e-10, info = paste("grad parity", bwtype))
  }
})
