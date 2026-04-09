test_that("mpi.cart.get parser validates output shape", {
  parse.fun <- getFromNamespace(".npRmpi_parse_cart_get", "npRmpi")

  parsed <- parse.fun(1:6, 2L)
  expect_identical(parsed$dims, 1:2)
  expect_identical(parsed$periods, 3:4)
  expect_identical(parsed$coords, 5:6)

  expect_error(parse.fun(1:5, 2L), "expected 6 values")
  expect_error(parse.fun(1:6, 0L), "single positive integer")
})

test_that("tail slave index helper stays in bounds for small sizes", {
  idx.fun <- getFromNamespace(".npRmpi_tail_slave_indices", "npRmpi")

  expect_identical(idx.fun(2L), 1L)
  expect_identical(idx.fun(3L), 1:2)
  expect_identical(idx.fun(9L), 7:8)
  expect_error(idx.fun(1L), "integer >= 2")
})
