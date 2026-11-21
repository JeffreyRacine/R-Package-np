test_that("degree-zero beta NN first derivatives match the frozen-radius reference", {
  if (!spawn_mpi_slaves(n = 1L)) skip("MPI workers unavailable")
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  # Saved minimum-count cases, with independent 80-digit derivative rows.
  cases <- readRDS(test_path("fixtures", "beta-nn-lc-derivative.rds"))
  for (v in cases) {
    x <- data.frame(x = v$x)
    grid <- data.frame(x = v$query)
    bw <- npregbw(xdat = x, ydat = v$y, bws = 2, bandwidth.compute = FALSE,
                  regtype = "lp", degree = 0L, bwtype = v$type, ckertype = "beta",
                  ckerbound = "fixed", ckerlb = v$bounds[1L], ckerub = v$bounds[2L])
    fit <- npreg(bw, txdat = x, tydat = v$y, exdat = grid, gradients = TRUE, se = FALSE)
    hat <- npreghat(bw, txdat = x, exdat = grid, s = 1L)
    expect_lt(abs(fit$grad[1L, 1L] - v$derivative), 1e-8)
    expect_identical(dim(hat), c(1L, length(v$x)))
    expect_lt(max(abs(as.numeric(hat) - v$hat)), 1e-8)
  }
})
