test_that("updateBwNameMetadata updates bandwidth metadata deterministically", {
  bws <- list(varnames = list(x = "oldx", y = "oldy"))
  names_in <- list(xnames = "x_new", ynames = "y_new")

  out <- updateBwNameMetadata(nameList = names_in, bws = bws)

  expect_equal(out$xnames, "x_new")
  expect_equal(out$ynames, "y_new")
  expect_equal(out$varnames$x, "x_new")
  expect_equal(out$varnames$y, "y_new")
  expect_equal(bws$varnames$x, "oldx")
  expect_equal(bws$varnames$y, "oldy")
})

test_that("updateBwNameMetadata is a no-op for empty metadata input", {
  bws <- list(varnames = list(x = "x0"))
  out <- updateBwNameMetadata(nameList = list(), bws = bws)
  expect_equal(out, bws)
})
