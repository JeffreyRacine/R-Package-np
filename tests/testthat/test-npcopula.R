test_that("npcopula basic functionality works", {
  data("faithful")
  # npcopula can take a dbandwidth (from npudistbw) for copula 
  # or bandwidth (from npudensbw) for copula density
  
  bw <- npudistbw(dat=faithful, bws=c(0.5, 5), bandwidth.compute=FALSE)
  
  # Copula CDF
  cop <- npcopula(data=faithful, bws=bw)
  expect_s3_class(cop, "data.frame")
  expect_true("copula" %in% names(cop))
  expect_true(all(cop$copula >= 0 & cop$copula <= 1))
  
  # Copula density
  bw_dens <- npudensbw(dat=faithful, bws=c(0.5, 5), bandwidth.compute=FALSE)
  cop_dens <- npcopula(data=faithful, bws=bw_dens)
  expect_s3_class(cop_dens, "data.frame")
  expect_true("copula" %in% names(cop_dens))
  expect_true(all(cop_dens$copula >= 0))
})
