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

# Evaluation admission is distinct from the radius/fit identity above.
test_that("fractional extended NN objectives share their integer cell", {
  withr::local_options(np.messages = FALSE, np.tree = FALSE, np.extendednn = TRUE)
  set.seed(2701)
  x <- data.frame(x = rnorm(60), z = rnorm(60))
  y <- sin(x$x) + .5*x$z + rnorm(60)
  for (type in c("generalized_nn", "adaptive_nn")) {
    evaluate <- function(k) {
      b <- npregbw(xdat = x, ydat = y, bws = c(k, 22), bwtype = type,
                   regtype = "lp", degree = c(2L, 2L), bandwidth.compute = FALSE)
      .npregbw_eval_only(x, y, b, invalid.penalty = "dbmax")$objective
    }
    for (k in c(59.5, 59.6, 60.4, 60.5, 60.6, 61.4)) {
      expected <- evaluate(round(k))
      expect_true(is.finite(expected) && abs(expected) < 1e100)
      expect_identical(evaluate(k), expected)
    }
  }
})
