test_that("npregiv basic functionality works", {
  set.seed(42)
  n <- 100
  # Simple IV setup
  z <- runif(n) # instrument
  w <- z + rnorm(n, sd=0.1) # endogenous regressor
  y <- w^2 + rnorm(n, sd=0.1) # outcome
  
  # This can be slow and sensitive to parameters
  # Use simple settings
  # Note: npregiv doesn't have a bandwidth.compute=FALSE in the same way, 
  # but we can try to keep it simple.
  
  # For the sake of a fast unit test, maybe we just check if it runs
  # without error with very small number of iterations if applicable,
  # but npregiv uses Tikhonov or Landweber-Fridman.
  
  # Let's try a very small example
  model <- npregiv(y=y, w=w, z=z, method="Landweber-Fridman", iterate.max=10)
  
  expect_type(model, "list")
  expect_true("phi" %in% names(model))
})

test_that("npregivderiv basic functionality works", {
  set.seed(42)
  n <- 50
  z <- runif(n)
  w <- z + rnorm(n, sd=0.1)
  y <- z^2 + rnorm(n, sd=0.1)
  
  # Just check if it runs
  model <- npregivderiv(y=y, z=z, w=w, iterate.max=2)
  expect_type(model, "list")
  expect_true("phi.prime" %in% names(model))
})
