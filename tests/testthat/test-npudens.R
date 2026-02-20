test_that("npudens basic functionality works", {
  # skip_on_cran()
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")

  data("faithful")
  # Use fixed bandwidths to avoid time-consuming cross-validation in tests
  bw <- npudensbw(dat=faithful, bws=c(0.5, 5), bandwidth.compute=FALSE)
  
  fit <- npudens(bws=bw)
  
  expect_s3_class(fit, "npdensity")
  expect_type(predict(fit), "double")
  expect_equal(length(predict(fit)), nrow(faithful))
  
  # Check summary and print
  expect_output(print(fit))
  expect_output(summary(fit))
})

test_that("npudens works with formula and factors", {
  # skip_on_cran()
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")

  data("Italy")
  Italy_sub <- Italy[1:100, ]
  bw <- npudensbw(formula=~ordered(year)+gdp, data=Italy_sub, 
                                bws=c(0.5, 1.0), bandwidth.compute=FALSE)
  
  fit <- npudens(bws=bw)
  expect_s3_class(fit, "npdensity")
})