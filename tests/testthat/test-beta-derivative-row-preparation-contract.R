locate_beta_derivative_row_sources <- function() {
  roots <- unique(c(
    Sys.getenv("NP_SOURCE_ROOT", unset = ""),
    normalizePath(file.path(getwd(), "..", ".."), mustWork = FALSE),
    normalizePath(getwd(), mustWork = FALSE)
  ))
  roots <- roots[nzchar(roots)]
  for (root in roots) {
    if (file.exists(file.path(root, "src", "continuous_kernel_row.c")) &&
        file.exists(file.path(root, "src", "continuous_kernel_row_gradient.c")))
      return(root)
  }
  NULL
}

test_that("beta derivative state shares the canonical PDF row preparation", {
  root <- locate_beta_derivative_row_sources()
  skip_if(is.null(root), "package sources unavailable")
  row <- paste(readLines(
    file.path(root, "src", "continuous_kernel_row.c"), warn = FALSE
  ), collapse = "\n")
  gradient <- paste(readLines(
    file.path(root, "src", "continuous_kernel_row_gradient.c"), warn = FALSE
  ), collapse = "\n")
  header <- paste(readLines(
    file.path(root, "src", "continuous_kernel_row.h"), warn = FALSE
  ), collapse = "\n")

  expect_false(grepl(
    "np_continuous_kernel_beta_prepared_derivative_row_prepare",
    paste(row, gradient, header), fixed = TRUE
  ))
  expect_equal(sum(gregexpr(
    "if(context->pdf_row_derivative_component != NULL) {",
    row, fixed = TRUE
  )[[1L]] > 0L), 1L)
  expect_match(
    gradient,
    paste0(
      "np_continuous_kernel_beta_prepared_derivative_context_prepare(\n",
      "        plan->beta_prepared);"
    ),
    fixed = TRUE
  )
  expect_match(
    gradient,
    paste0(
      "np_continuous_kernel_beta_prepared_pdf_row_prepare(\n",
      "      plan, segment, evaluation_index, omitted_observation,"
    ),
    fixed = TRUE
  )
})
