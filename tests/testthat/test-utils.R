test_that("npksum basic functionality works", {
  data("cps71")
  cps71_sub <- cps71[1:50,]
  
  k <- npksum(txdat=cps71_sub$age, tydat=cps71_sub$logwage, bws=1.0)
  
  expect_s3_class(k, "npkernelsum")
  expect_type(k$ksum, "double")
})

test_that("npseed works", {
  # npseed sets seed for C backend, doesn't affect R's runif
  expect_silent(npseed(42))
})

test_that("b.star works", {
  data("Italy")
  # b.star is for block size selection in bootstrap
  b <- b.star(Italy$gdp)
  expect_type(b, "double")
  expect_equal(nrow(b), 1)
})
