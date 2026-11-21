test_that("complete Gaussian extrapolation reproduces the fitted polynomial and its first derivative", {
  x <- data.frame(x = (-31:32) / 32)
  h <- 0.25
  bw <- npregbw(xdat = x, ydat = x$x, bws = h, bandwidth.compute = FALSE)
  pilot <- list(x = x, xcon = x, numeric.cols = 1L, factor.cols = integer(),
                groups = list(all = seq_len(nrow(x))), d = 1L, q = 5L,
                b = h, bernstein = TRUE, gaussian.pilot.range = range(x$x),
                pilot.bws = .np_bias_regbw(bw, x, h, 5L, TRUE))
  row <- .np_bias_pilot_rows(pilot, data.frame(x = -1.75), TRUE, 1L)
  # Independent 80-decimal weighted polynomial calculation on the binary inputs.
  donors <- c(1L, 2L, 4L, 8L, 16L, 32L)
  expected.mean <- c(3356.9387461917472, -4730.1948546369567, 1171.5314938711758, -300.67724462221798, 60.751883624178305, -0.07954286842600411)
  expected.derivative <- c(-23519.298351809706, 38628.669925944836, -13336.651739552944, 4077.8946383515527, -310.64148376437436, -0.15689906269113493)
  expect_lt(max(abs(row$mean[donors] - expected.mean)), 1e-8)
  expect_lt(max(abs(row$derivatives[[1L]][donors] - expected.derivative)), 1e-8)
})

test_that("Gaussian residual centering subtracts a polynomial on tied training coordinates", {
  x <- rep((-4:4) / 4, each = 2L)
  y <- 2^18 + sin(seq_along(x))
  r <- .np_bias_gaussian_center(x, y, 3L, range(x))
  removed <- y - r
  expect_identical(removed[seq(1,18,2)], removed[seq(2,18,2)])
  expect_lt(max(abs(diff(removed[seq(1,18,2)], differences = 4L))), 1e-8)
})

test_that("the Gaussian covariance correction uses the effective auxiliary degree", {
  x <- data.frame(x = seq(-0.97, 0.96, length.out = 120))
  y <- sin(3 * x$x) + 0.1 * cos(seq_along(x$x))
  bw <- npregbw(xdat = x, ydat = y, bws = 35L, bandwidth.compute = FALSE,
                regtype = "lc", bwtype = "generalized_nn", ckertype = "gaussian",
                ckerbound = "fixed", ckerlb = -1, ckerub = Inf)
  prepare <- function(control)
    .np_bias_regression_prepare(bw, x, y, control, adapt = FALSE)
  ordinary <- prepare(np_bias_control())
  explicit <- prepare(np_bias_control(degree = 5L))
  expect_identical(ordinary[["gaussian.pilot.range"]], range(ordinary$xcon[[1L]]))
  expect_identical(ordinary$residual, explicit$residual)
  expect_identical(ordinary[["gaussian.pilot.range"]], explicit[["gaussian.pilot.range"]])
  expect_null(prepare(np_bias_control(degree = 3L))[["gaussian.pilot.range"]])
})
