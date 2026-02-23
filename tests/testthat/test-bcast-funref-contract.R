test_that(".npRmpi_bcast_cmd_funref resolves standard command heads", {
  expect_identical(.npRmpi_bcast_cmd_funref(quote(assign)), "assign")
  expect_identical(.npRmpi_bcast_cmd_funref("assign"), "assign")
  expect_identical(.npRmpi_bcast_cmd_funref(base::assign), base::assign)
})

test_that(".npRmpi_bcast_cmd_funref resolves namespace-qualified heads", {
  ref <- .npRmpi_bcast_cmd_funref(quote(base::assign))
  expect_true(is.function(ref))
  expect_identical(ref, base::assign)

  ref_get <- .npRmpi_bcast_cmd_funref(quote(base::get))
  expect_true(is.function(ref_get))
  expect_identical(ref_get, base::get)
})
