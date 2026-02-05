test_that("npindex basic functionality works", {
  set.seed(42)
  n <- 100
  x1 <- runif(n)
  x2 <- runif(n)
  # Single index model: y = g(x1 + x2) + e
  y <- (x1 + x2)^2 + rnorm(n, sd=0.1)
  
  # Optimization might be slow, but for n=100 it should be okay.
  # Using a formula interface
  bw <- npindexbw(y~x1+x2, method="ichimura", nmulti=1)
  
  model <- npindex(bws=bw)
  
  expect_s3_class(model, "singleindex")
  expect_type(predict(model), "double")
  expect_output(summary(model))
})
