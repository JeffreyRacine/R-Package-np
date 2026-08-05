test_that("blocked bivariate entropy integration preserves the scalar formula", {
  helper <- getFromNamespace(".np_entropy_bivariate_integral", "npRmpi")
  set.seed(20260717)
  n <- 24L
  x <- rnorm(n)
  y <- 0.55 * x + sqrt(1 - 0.55^2) * rnorm(n)
  bw.x <- 0.42
  bw.y <- 0.47
  bw.joint <- c(0.50, 0.53)
  lower <- c(min(x) - 10 * IQR(x), min(y) - 10 * IQR(y))
  upper <- c(max(x) + 10 * IQR(x), max(y) + 10 * IQR(y))

  scalar.integrand <- function(xy) {
    f.x <- mean(dnorm((xy[1L] - x) / bw.x)) / bw.x
    f.y <- mean(dnorm((xy[2L] - y) / bw.y)) / bw.y
    f.xy <- mean(
      dnorm((xy[1L] - x) / bw.joint[1L]) *
        dnorm((xy[2L] - y) / bw.joint[2L]) /
        (bw.joint[1L] * bw.joint[2L])
    )
    (sqrt(f.xy) - sqrt(f.x) * sqrt(f.y))^2
  }
  reference <- 0.5 * cubature::adaptIntegrate(
    scalar.integrand, lowerLimit = lower, upperLimit = upper
  )$integral
  candidate <- helper(x, y, bw.x, bw.y, bw.joint, lower, upper)

  tolerance <- 64 * .Machine$double.eps * max(1, abs(reference))
  expect_equal(candidate, reference, tolerance = tolerance)
})

test_that("native entropy callback preserves pointwise Gaussian arithmetic", {
  set.seed(20260718)
  x <- rnorm(20)
  y <- rnorm(20)
  bandwidths <- c(0.4, 0.5, 0.45, 0.55)
  points <- cbind(
    c(0, 0),
    c(-1.25, 0.75),
    c(2.5, -3),
    c(-15, 15),
    c(25, -25)
  )

  scalar <- function(xy) {
    f.x <- mean(dnorm((xy[1L] - x) / bandwidths[1L])) / bandwidths[1L]
    f.y <- mean(dnorm((xy[2L] - y) / bandwidths[2L])) / bandwidths[2L]
    f.xy <- mean(
      dnorm((xy[1L] - x) / bandwidths[3L]) *
        dnorm((xy[2L] - y) / bandwidths[4L]) /
        (bandwidths[3L] * bandwidths[4L])
    )
    (sqrt(f.xy) - sqrt(f.x) * sqrt(f.y))^2
  }

  reference <- matrix(apply(points, 2L, scalar), nrow = 1L)
  candidate <- .Call(
    "C_np_entropy_gaussian_integrand",
    points,
    x,
    y,
    bandwidths,
    PACKAGE = "npRmpi"
  )

  tolerance <- 128 * .Machine$double.eps * max(1, abs(reference))
  expect_identical(dim(candidate), c(1L, ncol(points)))
  expect_equal(candidate, reference, tolerance = tolerance)
})

test_that("native entropy callback rejects malformed inputs without fallback", {
  call.native <- function(points, x = c(0, 1), y = c(0, 1),
                          bandwidths = rep(0.5, 4L)) {
    .Call(
      "C_np_entropy_gaussian_integrand",
      points,
      x,
      y,
      bandwidths,
      PACKAGE = "npRmpi"
    )
  }

  expect_error(call.native(c(0, 1)), "numeric matrix")
  expect_error(call.native(matrix(0, nrow = 3L)), "two rows")
  expect_error(call.native(matrix(0, nrow = 2L), y = 0), "equal positive")
  expect_error(
    call.native(matrix(0, nrow = 2L), bandwidths = c(0.5, 0, 0.5, 0.5)),
    "finite and positive"
  )
})
