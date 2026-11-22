test_that("linear beta ANN extrapolation recovers both reflected means and slopes", {
  old <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old), add = TRUE)
  # Physical-coordinate weighted fits computed at 80 digits with frozen donor radii.
  v <- readRDS(test_path("fixtures", "beta-ann-linear-upper.rds"))
  for (side in c(1, -1)) {
    x <- data.frame(x = side * v$x)
    grid <- data.frame(x = side * v$at)
    for (bernstein in c(FALSE, TRUE)) {
      bw <- npregbw(xdat = x, ydat = v$y, bws = v$count, bandwidth.compute = FALSE,
                    regtype = "lp", degree = 1L, basis = "glp", bernstein.basis = bernstein,
                    bwtype = "adaptive_nn", ckertype = "beta", ckerorder = 2L,
                    ckerbound = "fixed", ckerlb = min(side * v$bounds), ckerub = max(side * v$bounds))
      evaluate <- function() npreg(bw, txdat = x, tydat = v$y, exdat = grid,
                                   gradients = TRUE, se = FALSE)
      if (bernstein) {
        expect_warning(fit <- evaluate(), "outside training support")
      } else {
        fit <- evaluate()
      }
      expect_lt(max(abs(c(fit$mean, fit$grad[, 1L]) - c(1, side) * v$fit)), 1e-8)
      for (s in 0:1) {
        row <- npreghat(bw, txdat = x, exdat = grid, s = s, output = "matrix")
        expect_identical(dim(row), c(1L, length(v$x)))
        expect_equal(as.numeric(attr(row, "ridge.used", exact = TRUE)), 0)
        expect_lt(max(abs(as.numeric(row) - side^s * v$rows[s + 1L, ])), 1e-8)
      }
    }
  }
})
