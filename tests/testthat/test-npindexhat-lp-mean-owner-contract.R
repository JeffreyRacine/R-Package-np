library(np)

test_that("npindexhat lp mean owner preserves fit and matrix/apply parity for nonfixed bwtypes", {
  set.seed(20260315)
  n <- 90L
  x1 <- runif(n, -1, 1)
  x2 <- rnorm(n)
  y <- sin(2 * (0.7 * x1 - 0.3 * x2)) + 0.35 * x1 * x2 + rnorm(n, sd = 0.04)

  tx <- data.frame(x1 = x1, x2 = x2)
  ex <- data.frame(
    x1 = seq(min(x1) * 0.9, max(x1) * 0.9, length.out = 30L),
    x2 = seq(quantile(x2, 0.2), quantile(x2, 0.8), length.out = 30L)
  )

  cases <- list(
    list(label = "lp1 generalized canonical", bwtype = "generalized_nn", degree = 1L, basis = "glp", bern = FALSE, bws = c(1, 1, 9L)),
    list(label = "lp1 generalized tensor legacy", bwtype = "generalized_nn", degree = 1L, basis = "tensor", bern = FALSE, bws = c(1, 1, 9L)),
    list(label = "lp1 generalized bernstein legacy", bwtype = "generalized_nn", degree = 1L, basis = "glp", bern = TRUE, bws = c(1, 1, 9L)),
    list(label = "lp1 adaptive tensor", bwtype = "adaptive_nn", degree = 1L, basis = "tensor", bern = FALSE, bws = c(1, 1, 9L)),
    list(label = "lp2 generalized tensor", bwtype = "generalized_nn", degree = 2L, basis = "tensor", bern = FALSE, bws = c(1, 1, 11L)),
    list(label = "lp2 generalized bernstein", bwtype = "generalized_nn", degree = 2L, basis = "glp", bern = TRUE, bws = c(1, 1, 11L)),
    list(label = "lp2 adaptive tensor", bwtype = "adaptive_nn", degree = 2L, basis = "tensor", bern = FALSE, bws = c(1, 1, 11L))
  )

  for (cfg in cases) {
    bw <- npindexbw(
      xdat = tx,
      ydat = y,
      regtype = "lp",
      degree = cfg$degree,
      basis = cfg$basis,
      bernstein.basis = cfg$bern,
      bwtype = cfg$bwtype,
      bandwidth.compute = FALSE,
      bws = cfg$bws
    )

    fit.mean <- npindex(
      bws = bw,
      txdat = tx,
      tydat = y,
      exdat = ex,
      gradients = FALSE,
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

    expect_equal(as.vector(apply.mean), as.vector(fit.mean$mean), tolerance = 1e-8, info = paste("mean", cfg$label))
    expect_equal(as.vector(matrix.mean %*% y), as.vector(fit.mean$mean), tolerance = 1e-8, info = paste("mean matrix", cfg$label))
    expect_equal(as.vector(apply.mean), as.vector(matrix.mean %*% y), tolerance = 1e-10, info = paste("mean parity", cfg$label))
  }
})
