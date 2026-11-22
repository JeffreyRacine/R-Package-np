test_that("the linear beta ANN slope retains precision at the last training point", {
  opts <- options(np.messages = FALSE)
  on.exit(options(opts), add = TRUE)
  v <- readRDS(test_path("fixtures", "beta-ann-r112.rds"))
  x <- data.frame(x = v$x)
  bw <- npregbw(xdat = x, ydat = v$y, bws = v$count, bandwidth.compute = FALSE,
                regtype = "lp", degree = 1L, bernstein.basis = TRUE,
                ckertype = "beta", ckerorder = 2L, ckerbound = "fixed",
                ckerlb = v$bounds[1], ckerub = v$bounds[2], bwtype = "adaptive_nn")
  grid <- data.frame(x = max(v$x))
  H <- npreghat(bw, txdat = x, exdat = grid, s = 1L)
  expect_identical(dim(H), c(1L, nrow(x)))
  expect_equal(as.numeric(attr(H, "ridge.used", exact = TRUE)), 0)
  expect_lt(max(abs(as.numeric(H) - v$direction)), 1e-8)
  # A small row error is visible in slope variance even with unit donor variance.
  expect_lt(abs(sum(H * H) - v$unit_variance), 1e-8)
})
