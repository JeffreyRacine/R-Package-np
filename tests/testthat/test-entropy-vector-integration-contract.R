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

test_that("entropy integration block sizing is bounded and result invariant", {
  block.size <- getFromNamespace(".np_entropy_integration_block_size", "npRmpi")
  helper <- getFromNamespace(".np_entropy_bivariate_integral", "npRmpi")

  expect_identical(block.size(100L), 64L)
  expect_lt(block.size(100000L), 64L)
  expect_lte(64 * 100000 * block.size(100000L), 32 * 1024^2)

  set.seed(20260718)
  x <- rnorm(20)
  y <- rnorm(20)
  bw <- c(0.4, 0.5)
  lower <- c(min(x) - 10 * IQR(x), min(y) - 10 * IQR(y))
  upper <- c(max(x) + 10 * IQR(x), max(y) + 10 * IQR(y))
  ordinary <- helper(x, y, bw[1L], bw[2L], bw, lower, upper)
  forced <- helper(
    x, y, bw[1L], bw[2L], bw, lower, upper,
    target.bytes = 64 * length(x) * 4L
  )

  tolerance <- 64 * .Machine$double.eps * max(1, abs(ordinary))
  expect_equal(forced, ordinary, tolerance = tolerance)
})
