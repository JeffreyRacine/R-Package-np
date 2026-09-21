test_that("IV derivative initialization validates ties before preparation", {
  z <- c(0, 0, .3, .6, 1, 1)
  w <- seq_along(z)/length(z)
  y <- z^2+w/10
  sentinel <- function(...) stop("preparation reached")
  testthat::local_mocked_bindings(npudensbw=sentinel, .package="np")
  set.seed(4051)
  seed <- .Random.seed
  expect_error(npregivderiv(y=y,z=z,w=w,starting.values=seq_along(z)),
               "conflicting quadrature ordinates at duplicate coordinates")
  expect_identical(.Random.seed,seed)
  expect_error(npregivderiv(y=y,z=z,w=w,starting.values=2*z),
               "preparation reached")
  expect_error(npregivderiv(y=y,z=z,w=w), "preparation reached")
})
