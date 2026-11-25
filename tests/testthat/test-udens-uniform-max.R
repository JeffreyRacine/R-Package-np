test_that("the uniform Gaussian maximum respects shared innovations across targets", {
  F <- matrix(c(1, .7, -.2, .3, -.8, .6), nrow = 3)
  B <- 271L
  # The owner receives a square covariance factor, including any zero columns.
  factor <- cbind(F, 0)
  expected <- .np_with_seed(121L, {
    innovations <- matrix(rnorm(nrow(factor)*B), nrow(factor))
    apply(abs(factor %*% innovations), 2L, max)
  })
  actual <- .np_with_seed(121L, .np_udens_gaussian_absolute_maximum(factor, B, function(...) NULL))
  expect_equal(actual, expected, tolerance = 1e-14)
})
test_that("the joint envelope compares two marginal maximum quantiles", {
  a <- diag(c(1, .1))
  b <- diag(c(.1, 1))
  B <- 501L
  maxima <- .np_with_seed(27L, .np_udens_gaussian_absolute_maximum(a, B, function(...) NULL, b))
  z <- .np_with_seed(27L, matrix(rnorm(2L * B), 2L))
  reference <- cbind(apply(abs(a %*% z), 2L, max), apply(abs(b %*% z), 2L, max))
  expect_equal(unname(maxima), reference, tolerance = 1e-14)
  expect_lt(max(apply(maxima, 2L, quantile, .95)), unname(quantile(apply(maxima, 1L, max), .95)))
  g <- rep(seq(-1, 1, length.out = 11), each = 2)
  factor <- .np_udens_uniform_empirical_factor(seq(-2, 2, length.out = 9), .43, g)
  expect_equal(factor, t(factor), tolerance = 1e-12)
  expect_equal(factor[1, ], factor[2, ], tolerance = 1e-12)
})

test_that("uniform joint critical value and multiblock influence include the fitted guide", {
  x <- seq(-2.2, 2.5, length.out = 83L) + .09 * sin(seq_len(83L))
  g <- c(-1.1, -.7, -.15, .32, .74)
  h <- .52
  rule <- .np_udens_q6_rule(64L)
  kernel <- function(u) ifelse(abs(u) < 1, (225 - 1050*u^2 + 945*u^4)/128, 0)
  k <- kernel(outer(x, g, "-") / h) / h
  nodes <- outer(g, h * rule$x, "+")
  weights <- kernel(rule$x) * rule$w
  functional <- function(p) {
    location <- sum(p*x)
    sigma <- sqrt(sum(p*(x-location)^2))
    drop(crossprod(p,k)) + dnorm(g,location,sigma) -
      drop(dnorm(nodes,location,sigma) %*% weights)
  }
  n <- length(x)
  p <- rep(1/n,n)
  eps <- 2e-6
  rows <- t(vapply(seq_len(n), function(i) {
    towards <- -p
    towards[i] <- towards[i] + 1
    (functional(p+eps*towards) - functional(p-eps*towards))/(2*eps)
  }, numeric(length(g))))
  F <- .np_udens_uniform_empirical_factor(x,h,g)
  expect_lt(max(abs(tcrossprod(F)-crossprod(rows)/n^2)),1e-10)
  state <- list(construction="density-uniform-q6-normal-v1",x=data.frame(x=x),
                n=n,h=h,plot.random.seed=309L)
  band <- .np_udens_uniform_normal_public(state,data.frame(x=g),.05,"all",9999L)$inference
  critical <- max(band$settings$critical.components)
  expect_equal(band$critical,critical,tolerance=1e-14)
  expect_equal(as.numeric(band$bands$simultaneous[,2L,1L]),band$center+critical,tolerance=1e-14)
})
