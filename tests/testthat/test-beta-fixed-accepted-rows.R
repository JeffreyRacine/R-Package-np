test_that("fixed beta fits and hats recover the accepted local polynomial near a boundary", {
  opts <- options(np.messages = FALSE)
  on.exit(options(opts), add = TRUE)
  v <- readRDS(test_path("fixtures", "beta-fixed-r112.rds"))
  x <- data.frame(x = v$x)
  for (m in v$models) {
    bw <- npregbw(xdat = x, ydat = v$y, bws = m$bw, bandwidth.compute = FALSE,
                  regtype = "lp", degree = m$degree, bernstein.basis = TRUE,
                  ckertype = "beta", ckerorder = 2L, ckerbound = "fixed",
                  ckerlb = m$bounds[1], ckerub = m$bounds[2])
    grid <- data.frame(x = m$grid)
    fit <- npreg(bw, txdat = x, tydat = v$y, exdat = grid, gradients = TRUE, se = FALSE)
    for (s in 0:1) {
      H <- npreghat(bw, txdat = x, exdat = grid, s = s)
      ref <- m$rows[[s + 1L]]
      expect_identical(dim(H), dim(ref))
      expect_identical(as.double(attr(H, "ridge.used", exact = TRUE)), numeric(nrow(grid)))
      expect_lt(max(abs(H - ref)), 1e-8)
      estimate <- if (s == 0L) fit$mean else fit$grad[, 1L]
      expect_lt(max(abs(estimate - drop(ref %*% v$y))), 1e-8)
    }
  }
})
