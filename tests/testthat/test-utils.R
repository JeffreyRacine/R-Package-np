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

test_that("b.star round=TRUE applies elementwise bounds for multivariate input", {
  set.seed(1)
  x <- cbind(arima.sim(n = 200, list(ar = 0.5)),
             arima.sim(n = 200, list(ar = 0.8)))

  b.raw <- b.star(x, round = FALSE)
  b.round <- b.star(x, round = TRUE)
  b.max <- ceiling(min(3 * sqrt(nrow(x)), nrow(x) / 3))

  expected <- cbind(
    pmax(1, pmin(b.max, round(b.raw[, "BstarSB"]))),
    pmax(1, pmin(b.max, round(b.raw[, "BstarCB"])))
  )
  colnames(expected) <- colnames(b.raw)

  expect_equal(b.round, expected)
})
