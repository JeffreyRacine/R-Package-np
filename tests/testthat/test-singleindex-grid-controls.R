test_that("single-index plots honor index quantiles and refuse slice holding", {
  d <- data.frame(x = seq(-1, 1, length.out = 30), z = cos(seq_len(30)))
  d$y <- sin(d$x) + .2 * d$z
  b <- npindexbw(y ~ x + z, data = d, bws = c(1, .5, .3), bandwidth.compute = FALSE)
  ix <- drop(as.matrix(d[c("x", "z")]) %*% b$beta)
  out <- plot(b, output = "data", neval = 7L,
              grid_control = np_grid_control(xtrim = c(.1, .8)))
  oracle <- seq(quantile(ix, .1), quantile(ix, .8), length.out = 7L)
  expect_equal(out$panels[[1L]]$grid[[1L]], oracle, tolerance = 1e-12)
  expect_error(plot(b, output = "data", grid_control = np_grid_control(xq = .2)),
               "xq is not applicable", fixed = TRUE)
})
