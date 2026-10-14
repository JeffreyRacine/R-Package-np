test_that("beta general-LP derivative HC0 covers every order and bandwidth mode", {
  old <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old), add = TRUE)

  n <- 15L
  txdat <- data.frame(
    x1 = seq(0.04, 0.96, length.out = n),
    x2 = 0.5 + 0.43 * sin(seq(0.2, 2.8, length.out = n))
  )
  tydat <- sin(4 * txdat$x1) + 0.35 * txdat$x2 + seq_len(n) / 140
  exdat <- txdat[c(2L, 5L, 8L, 12L, 15L), , drop = FALSE]

  for (order in c(2L, 4L, 6L, 8L)) {
    for (bwtype in c("fixed", "generalized_nn", "adaptive_nn")) {
      cell.x <- txdat
      cell.y <- tydat
      cell.ex <- exdat
      if (order == 2L && identical(bwtype, "adaptive_nn")) {
        # The original curved design is nearly interpolating at its last
        # row (q about 2.7e-13). Public unit responses change beta's retained
        # moment scale/forward rounding before q normalization, so that map
        # is not a portable strict numerical oracle for the original stress
        # fixture below. Retain the same kernel/order/bandwidth/basis contract
        # on a non-nearly-interpolating design for this reference comparison.
        cell.x$x2 <- (((seq_len(n) * 7L) %% n) + 0.5) / n
        cell.y <- sin(4 * cell.x$x1) + 0.35 * cell.x$x2 + seq_len(n) / 140
        cell.ex <- cell.x[c(2L, 5L, 8L, 12L, 15L), , drop = FALSE]
      }
      bw <- h6_explicit_lp_bw(
        cell.x, cell.y,
        if (identical(bwtype, "fixed")) c(0.27, 0.3) else c(7, 7),
        bwtype = bwtype, degree = c(2L, 1L),
        basis = "glp", bernstein = order %in% c(4L, 8L),
        ckertype = "beta", ckerorder = order,
        ckerbound = "fixed", ckerlb = c(0, 0), ckerub = c(1, 1)
      )
      h6_expect_lp_derivative_contract(
        bw, cell.x, cell.y, cell.ex, tolerance = 1e-5,
        oracle.method = "actual"
      )
    }
  }
})

