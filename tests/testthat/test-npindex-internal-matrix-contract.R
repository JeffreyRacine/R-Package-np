test_that("npindex internal geometry consumes the canonical numeric matrix", {
  validate <- getFromNamespace(".npindex_internal_numeric_matrix", "npRmpi")
  index_from_beta <- getFromNamespace(".npindex_index_from_beta_tail", "npRmpi")
  coordinate_setup <- getFromNamespace(".npindex_beta_coordinate_setup", "npRmpi")
  xmat <- matrix(seq_len(18L) / 10, ncol = 3L)

  expect_identical(validate(xmat, "test"), xmat)
  expect_error(
    validate(as.data.frame(xmat), "test"),
    "requires a numeric matrix",
    fixed = TRUE
  )
  expect_error(
    validate(matrix(letters[seq_along(xmat)], ncol = 3L), "test"),
    "requires a numeric matrix",
    fixed = TRUE
  )

  index.source <- paste(deparse(index_from_beta), collapse = "\n")
  coordinate.source <- paste(deparse(coordinate_setup), collapse = "\n")
  expect_match(index.source, ".npindex_internal_numeric_matrix", fixed = TRUE)
  expect_match(coordinate.source, ".npindex_internal_numeric_matrix", fixed = TRUE)
  expect_false(grepl("toMatrix", index.source, fixed = TRUE))
  expect_false(grepl("toMatrix", coordinate.source, fixed = TRUE))
})
