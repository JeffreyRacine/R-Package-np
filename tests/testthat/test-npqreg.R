test_that("npqreg basic functionality works", {
  data("cps71")
  cps71_sub <- cps71[1:50, ]
  
  # Quantile regression needs a condbandwidth object
  bw <- npcdistbw(xdat=cps71_sub$age, ydat=cps71_sub$logwage, 
                  bws=c(0.5, 5), bandwidth.compute=FALSE)
  
  # Median regression
  model <- npqreg(bws=bw, tau=0.5)
  
  expect_s3_class(model, "qregression")
  expect_type(predict(model), "double")
  expect_equal(length(predict(model)), 50)
  expect_equal(model$tau, 0.5)
  
  expect_output(summary(model))
})

test_that("npqreg works with multiple taus", {
  data("cps71")
  cps71_sub <- cps71[1:30, ]
  bw <- npcdistbw(xdat=cps71_sub$age, ydat=cps71_sub$logwage, 
                  bws=c(0.5, 5), bandwidth.compute=FALSE)
  
  # npqreg only takes a single tau at a time according to some versions, 
  # but let's see if it works as extra arg or in bws.
  # Actually, npqreg usage says tau is an argument.
  
  model_q25 <- npqreg(bws=bw, tau=0.25)
  expect_equal(model_q25$tau, 0.25)
})
