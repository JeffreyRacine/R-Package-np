test_that("npksum basic functionality works", {
  # # skip_on_cran()
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  
  data("cps71")
  cps71_sub <- cps71[1:50,]
  k <- npksum(txdat=cps71_sub$age, tydat=cps71_sub$logwage, bws=1.0)
  
  expect_s3_class(k, "npkernelsum")
  expect_type(k$ksum, "double")
})

test_that("npseed works", {
  # # skip_on_cran()
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  
  # npseed sets seed for C backend
  expect_silent(npseed(42))
})

test_that("b.star works", {
  # b.star is a standard R function in npRmpi, doesn't necessarily need MPI 
  # but we can run it on master.
  data("Italy")
  b <- b.star(Italy$gdp)
  expect_type(b, "double")
  expect_equal(nrow(b), 1)
})