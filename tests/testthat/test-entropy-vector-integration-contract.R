test_that("deterministic bivariate entropy quadrature preserves its oracle", {
  helper <- getFromNamespace(".np_entropy_bivariate_integral", "np")
  domain.helper <- getFromNamespace(".np_entropy_bivariate_domain", "np")
  set.seed(20260717)
  n <- 24L
  x <- rnorm(n)
  y <- 0.55 * x + sqrt(1 - 0.55^2) * rnorm(n)
  bw.x <- 0.42
  bw.y <- 0.47
  bw.joint <- c(0.50, 0.53)
  domain <- domain.helper(x, y, bw.x, bw.y, bw.joint)
  candidate <- helper(
    x, y, bw.x, bw.y, bw.joint, domain$lower, domain$upper
  )

  ## Independently evaluated at tighter two-dimensional quadrature settings.
  reference <- 0.008252664978
  expect_lte(abs(candidate - reference), 1e-8)
})

test_that("entropy trapezoid weights match the established helper", {
  weights <- getFromNamespace(".np_entropy_trapezoid_weights", "np")
  established <- getFromNamespace("integrate.trapezoidal", "np")
  grid <- seq(-4, 5, length.out = 257L)
  value <- exp(-grid^2 / 3) * (1 + cos(grid)^2)

  expect_equal(
    unname(crossprod(weights(grid), value)[1L]),
    established(grid, value)[length(grid)],
    tolerance = 32 * .Machine$double.eps
  )
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
    PACKAGE = "np"
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
      PACKAGE = "np"
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

test_that("entropy acceleration preserves the explicit scalar-option contract", {
  set.seed(20260808)
  x <- rnorm(64)
  y <- rnorm(64)
  points <- rbind(seq(-2, 2, length.out = 31L),
                  seq(2, -2, length.out = 31L))
  bandwidths <- c(0.4, 0.5, 0.45, 0.55)
  old <- options(np.macMseries.accelerate = FALSE)
  on.exit(options(old), add = TRUE)

  scalar <- .Call("C_np_entropy_gaussian_integrand", points, x, y,
                  bandwidths, PACKAGE = "np")
  options(np.macMseries.accelerate = TRUE)
  accelerated <- .Call("C_np_entropy_gaussian_integrand", points, x, y,
                       bandwidths, PACKAGE = "np")
  options(np.macMseries.accelerate = FALSE)
  scalar.repeat <- .Call("C_np_entropy_gaussian_integrand", points, x, y,
                         bandwidths, PACKAGE = "np")

  expect_identical(scalar.repeat, scalar)
  expect_equal(accelerated, scalar,
               tolerance = 128 * .Machine$double.eps)
})
