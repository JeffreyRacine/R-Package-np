test_that("npcmstest basic functionality works", {
  set.seed(42)
  n <- 50
  x <- rnorm(n)
  y <- 1 + x + rnorm(n, sd=0.1)
  # Parametric model - need x=TRUE, y=TRUE for npcmstest
  model <- lm(y~x, x=TRUE, y=TRUE)
  
  # Asymptotic for speed
  test <- npcmstest(model=model, xdat=x, ydat=y, distribution="asymptotic")
  
  expect_s3_class(test, "cmstest")
  expect_output(summary(test))
})

test_that("npqcmstest basic functionality works", {
  library(quantreg)
  set.seed(42)
  n <- 50
  x <- rnorm(n)
  y <- 1 + x + rnorm(n, sd=0.1)
  # npqcmstest needs an rq model with model=TRUE
  model <- rq(y~x, tau=0.5, model=TRUE)
  
  # Asymptotic for speed
  test <- npqcmstest(model=model, xdat=x, ydat=y, distribution="asymptotic")
  
  expect_s3_class(test, "cmstest") 
})

test_that("npdeneqtest basic functionality works", {
  set.seed(42)
  n <- 50
  x <- data.frame(v1=rnorm(n))
  y <- data.frame(v1=rnorm(n, mean=0.5))
  
  # Use small boot.num
  test <- npdeneqtest(x, y, boot.num=19)
  
  expect_s3_class(test, "deneqtest")
  expect_output(summary(test))
})

test_that("npsymtest basic functionality works", {
  set.seed(42)
  n <- 50
  x <- rgamma(n, shape=2)
  
  test <- npsymtest(x, method="summation", boot.num=19)
  
  expect_s3_class(test, "symtest")
  expect_output(summary(test))
})

test_that("npunitest basic functionality works", {
  set.seed(42)
  n <- 50
  x <- rnorm(n)
  y <- rnorm(n, mean=0.5)
  
  test <- npunitest(x, y, method="summation", boot.num=19)
  
  expect_s3_class(test, "unitest")
  expect_output(summary(test))
})
