test_that("bounded adaptive quadratic hats and fitted derivatives agree with the saved weighted fits", {
  old <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old), add = TRUE)
  # Independent 80-digit local quadratic reference with the saved donor radii.
  cases <- readRDS(test_path("fixtures", "adaptive-quadratic-precision.rds"))
  for (v in cases) {
    x <- data.frame(x = v$x)
    grid <- data.frame(x = v$query)
    bw <- npregbw(xdat = x, ydat = v$y, bws = v$count, bandwidth.compute = FALSE,
                  regtype = "lp", degree = 2L, basis = "glp", bernstein.basis = v$bernstein,
                  bwtype = "adaptive_nn", ckertype = v$kernel, ckerorder = 2L,
                  ckerbound = "fixed", ckerlb = v$bounds[1L], ckerub = v$bounds[2L])
    fit <- npreg(bw, txdat = x, tydat = v$y, exdat = grid, gradients = TRUE, se = FALSE)
    expect_lt(max(abs(c(fit$mean, fit$grad[, 1L]) - v$fit)), 1e-8)
    for (s in 0:1) {
      row <- npreghat(bw, txdat = x, exdat = grid, s = s, output = "matrix")
      expect_equal(dim(row), c(1L, length(v$x)))
      expect_lt(max(abs(as.numeric(row) - v$rows[[s + 1L]])), 1e-8)
    }
  }
})
