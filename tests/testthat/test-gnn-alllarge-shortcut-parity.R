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


test_that("an aborted polynomial shortcut leaves no partial CV score", {
  old <- options(np.messages = FALSE, np.tree = FALSE,
                 np.extendednn = TRUE, np.macMseries.accelerate = FALSE)
  on.exit(options(old), add = TRUE)
  x <- c(rep(0, 25L), rep(1, 24L), 3)
  y <- x + sin(seq_along(x))
  .npRmpi_with_local_regression({
    for (position in c(1L, 25L, 50L)) {
      order <- append(seq_len(49L), 50L, after = position - 1L)
      data <- data.frame(x = x[order])
      response <- y[order]
      for (type in c("fixed", "generalized_nn")) {
        bw <- npregbw(xdat = data, ydat = response,
                      bws = if (type == "fixed") 100 else 5000,
                      bandwidth.compute = FALSE, bwtype = type,
                      regtype = "lp", degree = 2L,
                      bernstein.basis = FALSE)
        options(np.largeh = FALSE)
        literal <- .npregbw_eval_only(data, response, bw,
                                      invalid.penalty = "dbmax")
        options(np.largeh = TRUE)
        attempted <- .npregbw_eval_only(data, response, bw,
                                        invalid.penalty = "dbmax")
        expect_identical(as.integer(attempted$num.feval.fast), 0L)
        expect_equal(attempted$objective, literal$objective, tolerance = 1e-10)
      }
    }
  })
})
