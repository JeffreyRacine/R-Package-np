test_that("extended generalized-NN regression fits preserve training identity", {
  old <- options(np.messages = FALSE, np.tree = FALSE,
                 np.extendednn = TRUE)
  on.exit(options(old), add = TRUE)

  x <- data.frame(x = c(-1.8, -0.9, -0.25, 0.1, 0.65, 1.4, 2.3))
  y <- sin(1.3 * x$x) + seq_len(nrow(x)) / 50
  k <- nrow(x) + 2L
  bw <- npregbw(
    xdat = x, ydat = y, bws = k, bwmethod = "cv.ls",
    bwtype = "generalized_nn", bwscaling = FALSE,
    regtype = "lc", ckertype = "gaussian",
    bandwidth.compute = FALSE)

  expected <- vapply(seq_len(nrow(x)), function(index) {
    distance <- sort(abs(x$x[-index] - x$x[[index]]), method = "radix")
    radius <- max(distance) * k / length(distance)
    weight <- stats::dnorm((x$x[[index]] - x$x) / radius) / radius
    sum(weight * y) / sum(weight)
  }, numeric(1L))

  for (tree in c(FALSE, TRUE)) {
    options(np.tree = tree)
    expect_equal(
      fitted(npreg(bws = bw, txdat = x, tydat = y)),
      expected,
      tolerance = 2e-12)
  }
})
