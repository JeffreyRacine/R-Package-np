test_that("npcdist basic functionality works", {
  data("faithful")
  bw <- npcdistbw(xdat=faithful$eruptions, ydat=faithful$waiting, 
                  bws=c(0.5, 5), bandwidth.compute=FALSE)
  
  fit <- npcdist(bws=bw)
  
  expect_s3_class(fit, "condistribution")
  expect_type(predict(fit), "double")
  expect_equal(length(predict(fit)), nrow(faithful))
  
  expect_output(summary(fit))
})
