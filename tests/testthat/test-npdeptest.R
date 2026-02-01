test_that("npdeptest basic functionality works", {
  set.seed(42)
  n <- 50
  x <- rnorm(n)
  y <- x + rnorm(n, sd=0.1)
  
  # Use summation and small boot.num for speed
  test <- npdeptest(x, y, method="summation", boot.num=19)
  
  expect_s3_class(test, "deptest")
  expect_output(summary(test))
})

test_that("npsdeptest basic functionality works", {
  set.seed(42)
  n <- 50
  y <- arima.sim(n=n, list(ar=0.5))
  
  test <- npsdeptest(y, lag.num=1, method="summation", boot.num=19)
  
  expect_s3_class(test, "sdeptest")
  expect_output(summary(test))
})
