test_that(".force.type rejects unsupported type codes", {
  fn <- getFromNamespace(".force.type", "npRmpi")

  expect_identical(fn(1, 0), 1)
  expect_identical(fn(1, 6), 1)
  expect_error(fn(1, NA_integer_), "single non-missing integer code")
})

test_that("raw chunk splitter preserves zero-length chunks safely", {
  split.fun <- getFromNamespace(".npRmpi_split_raw_by_counts", "npRmpi")

  chunks <- split.fun(as.raw(1:3), c(1L, 0L, 2L), "unit-test")

  expect_type(chunks, "list")
  expect_identical(length(chunks), 3L)
  expect_identical(chunks[[1L]], as.raw(1))
  expect_identical(chunks[[2L]], raw(0))
  expect_identical(chunks[[3L]], as.raw(2:3))
})

test_that("raw chunk splitter rejects mismatched receive counts", {
  split.fun <- getFromNamespace(".npRmpi_split_raw_by_counts", "npRmpi")

  expect_error(split.fun(as.raw(1:2), c(1L, 2L), "unit-test"),
               "receive counts sum 3 but payload length is 2")
  expect_error(split.fun(as.raw(1:2), c(1L, NA_integer_), "unit-test"),
               "receive counts contain NA")
  expect_error(split.fun(as.raw(1:2), c(1L, -1L), "unit-test"),
               "receive counts must be non-negative")
})

test_that("raw length validator rejects invalid lengths", {
  validate.fun <- getFromNamespace(".npRmpi_validate_raw_length", "npRmpi")

  expect_identical(validate.fun(3L, "unit-test"), 3L)
  expect_error(validate.fun(NA_integer_, "unit-test"),
               "single non-missing length value")
  expect_error(validate.fun(-1L, "unit-test"),
               "non-negative length value")
})
