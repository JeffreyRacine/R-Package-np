test_that("npplreg basic functionality works", {
  set.seed(42)
  n <- 100
  x1 <- runif(n) # nonparametric part
  z1 <- runif(n) # parametric part
  y <- x1^2 + 2*z1 + rnorm(n, sd=0.1)
  
  # Partially linear model: y = g(x1) + z1*beta + e
  # bws needs to be a matrix. Row 1: y on x1, Row 2: z1 on x1
  bw_mat <- matrix(c(0.1, 0.1), nrow=2, ncol=1)
  bw <- npplregbw(xdat=z1, zdat=x1, ydat=y, bws=bw_mat, bandwidth.compute=FALSE)
  
  model <- npplreg(bws=bw)
  
  expect_s3_class(model, "plregression")
  expect_type(predict(model), "double")
  expect_output(summary(model))
  expect_length(coef(model), 1)
})
