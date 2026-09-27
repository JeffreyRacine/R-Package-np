test_that("direct regression fits the retained physical bandwidth", {

  old <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old), add = TRUE)
  set.seed(2201)
  x <- data.frame(x = rnorm(36), u = factor(rep(c("a", "b", "c"), 12)))
  y <- sin(x$x) + .3 * (x$u == "b")
  b <- npregbw(xdat = x, ydat = y, bws = c(1.2, .2),
               bwscaling = TRUE, bandwidth.compute = FALSE, regtype = "lc")
  h <- unlist(b$bandwidth)
  K <- dnorm(outer(x$x, x$x, "-") / h[1]) *
    ifelse(outer(x$u, x$u, "=="), 1 - h[2], h[2] / 2)
  expected <- colSums(K * y) / colSums(K)
  direct <- getFromNamespace(".np_regression_direct", "np")
  actual <- direct(b, x, y)
  expect_equal(actual$mean, expected, tolerance = 2e-12)
})
