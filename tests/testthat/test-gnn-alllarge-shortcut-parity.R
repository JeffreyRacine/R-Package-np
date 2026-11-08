test_that("GNN polynomial CV uses the certified global shortcut", {
  old <- options(np.messages = FALSE, np.tree = FALSE,
                 np.largeh = TRUE, np.extendednn = TRUE,
                 np.macMseries.accelerate = FALSE)
  on.exit(options(old), add = TRUE)
  x <- data.frame(x = seq(-1, 1, length.out = 48L))
  y <- sin(3 * x$x) + cos(7 * x$x) / 4
  for (degree in 1:3) {
    .npRmpi_with_local_regression({
      bw <- npregbw(xdat = x, ydat = y, bws = 4000,
                    bandwidth.compute = FALSE, bwtype = "generalized_nn",
                    regtype = "lp", degree = degree,
                    bernstein.basis = TRUE)
      options(np.largeh = TRUE)
      fast <- .npregbw_eval_only(x, y, bw, invalid.penalty = "dbmax")
      # Independent polynomial QR/hat-diagonal LOO, without package kernels.
      design <- outer(x$x, 0:degree, `^`)
      decomposition <- qr(design)
      estimate <- drop(design %*% qr.coef(decomposition, y))
      leverage <- rowSums(qr.Q(decomposition)^2)
      expected <- mean(((y - estimate) / (1 - leverage))^2)
      expect_equal(fast$objective, expected, tolerance = 1e-10)
      expect_identical(as.integer(fast$num.feval.fast), 1L)
      options(np.largeh = FALSE)
      literal <- .npregbw_eval_only(x, y, bw, invalid.penalty = "dbmax")
      expect_true(is.finite(literal$objective))
      expect_identical(as.integer(literal$num.feval.fast), 0L)
    })
  }
})
