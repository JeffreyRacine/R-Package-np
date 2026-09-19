test_that("density-equality observed activity preserves the public statistic", {
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  set.seed(193)
  x <- data.frame(x = rnorm(43))
  y <- data.frame(x = rnorm(27))
  result <- npdeneqtest(x, y, bw.x = .6, bw.y = .7, B = 9)
  sums <- function(a, b = a, h, loo = FALSE, power = 1L) {
    args <- list(txdat = a, bws = h, leave.one.out = loo,
                 bandwidth.divide = TRUE, kernel.pow = power)
    if (!loo) args$exdat <- b
    sum(do.call(npksum, args)$ksum)
  }
  n1 <- nrow(x); n2 <- nrow(y)
  In <- sums(x, h = .6, loo = TRUE) / (n1 * (n1 - 1)) +
    sums(y, h = .7, loo = TRUE) / (n2 * (n2 - 1)) -
    2 * sums(x, y, h = .6) / (n1 * n2)
  variance <- 2 * (sums(x, h = .6, loo = TRUE, power = 2L) / (n1^2 * (n1-1)^2) +
    sums(y, h = .7, loo = TRUE, power = 2L) / (n2^2 * (n2-1)^2) +
    2 * sums(x, y, h = .6, power = 2L) / (n1^2 * n2^2))
  expect_equal(result$In, In, tolerance = 2e-12)
  expect_equal(result$Tn, In / sqrt(variance), tolerance = 2e-12)
  expect_null(getFromNamespace(".np_progress_runtime", "np")$fit_state)
  expect_null(getFromNamespace(".np_progress_runtime", "np")$fit_forward)
})
