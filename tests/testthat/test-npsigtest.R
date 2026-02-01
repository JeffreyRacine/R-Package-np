test_that("npsigtest basic functionality works", {
  set.seed(42)
  n <- 50 # Keep it small for speed
  x1 <- runif(n)
  x2 <- runif(n)
  y <- x1^2 + rnorm(n, sd=0.1) # x2 is irrelevant
  
  bw <- npregbw(y~x1+x2, bws=c(0.1, 0.5), bandwidth.compute=FALSE)
  
  # Significance test can be slow, use few boot replications
  sig <- npsigtest(bws=bw, boot.num=19)
  
  expect_s3_class(sig, "sigtest")
  expect_output(summary(sig))
})
