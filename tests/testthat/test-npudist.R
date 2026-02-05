test_that("npudist basic functionality works", {
  data("faithful")
  # Use fixed bandwidths for speed
  bw <- npudistbw(dat=faithful, bws=c(0.5, 5), bandwidth.compute=FALSE)
  
  fit <- npudist(bws=bw)
  
  expect_s3_class(fit, "npdistribution")
  expect_type(predict(fit), "double")
  expect_true(all(predict(fit) >= 0 & predict(fit) <= 1))
  expect_equal(length(predict(fit)), nrow(faithful))
  
  expect_output(summary(fit))
})

test_that("npudist works with formula and factors", {
  data("Italy")
  Italy_sub <- Italy[1:100, ]
  bw <- npudistbw(formula=~ordered(year)+gdp, data=Italy_sub, 
                  bws=c(0.5, 1.0), bandwidth.compute=FALSE)
  
  fit <- npudist(bws=bw)
  expect_s3_class(fit, "npdistribution")
})
