test_that("closed npksum value calls reject lexical or ambiguous arguments", {
  closed_call <- getFromNamespace(".npksum_closed_values", "npRmpi")

  expect_error(
    closed_call(list(1)),
    "uniquely named list",
    fixed = TRUE
  )
  expect_error(
    closed_call(list(bws = 1, bws = 2)),
    "uniquely named list",
    fixed = TRUE
  )
  expect_error(
    closed_call(list(bws = quote(local_symbol))),
    "values, not language objects",
    fixed = TRUE
  )
})

test_that("reghat npksum builders use only the shared closed-value owner", {
  closed_call <- getFromNamespace(".npksum_closed_values", "npRmpi")
  owner <- getFromNamespace(
    ".npreghat_exact_lc_derivative_matrix_from_npksum_chunked",
    "npRmpi"
  )
  closed_text <- paste(
    deparse(body(closed_call), width.cutoff = 500L),
    collapse = "\n"
  )
  body_text <- paste(deparse(body(owner), width.cutoff = 500L), collapse = "\n")

  expect_true(grepl('as.name("npksum.default")', closed_text, fixed = TRUE))
  expect_false(grepl("quote(npksum.default", body_text, fixed = TRUE))
  expect_true(grepl(".npksum_closed_values", body_text, fixed = TRUE))
})
