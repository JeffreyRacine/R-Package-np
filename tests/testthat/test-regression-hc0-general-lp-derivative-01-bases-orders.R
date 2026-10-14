test_that("all LP bases and requested derivative orders share the HC0 adjoint block", {
  old <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old), add = TRUE)

  n <- 29L
  txdat <- data.frame(
    x1 = seq(0.03, 0.97, length.out = n),
    x2 = (((seq_len(n) * 7L) %% n) + 0.5) / n
  )
  tydat <- sin(2 * pi * txdat$x1) + 0.3 * txdat$x2^2 +
    seq_len(n) / 150
  exdat <- txdat[c(3L, 10L, 18L, 27L), , drop = FALSE]
  specifications <- list(
    list(basis = "glp", bernstein = FALSE),
    list(basis = "glp", bernstein = TRUE),
    list(basis = "tensor", bernstein = TRUE)
  )

  for (specification in specifications) {
    bw <- h6_explicit_lp_bw(
      txdat, tydat, c(0.3, 0.33), degree = c(2L, 2L),
      basis = specification$basis,
      bernstein = specification$bernstein,
      ckerbound = "fixed", ckerlb = c(0, 0), ckerub = c(1, 1)
    )
    h6_expect_lp_derivative_contract(
      bw, txdat, tydat, exdat, gradient.order = c(2L, 1L),
      tolerance = 2e-8, oracle.method = "actual"
    )
  }
})

