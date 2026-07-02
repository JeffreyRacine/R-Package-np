read_header_define <- function(symbol) {
  header <- testthat::test_path("..", "..", "src", "headers.h")
  skip_if_not(file.exists(header), "source headers.h is not available in this test context")
  lines <- readLines(header, warn = FALSE)
  pattern <- sprintf("^#define[[:space:]]+%s[[:space:]]+([0-9]+)\\b", symbol)
  hit <- grep(pattern, lines, value = TRUE)
  expect_length(hit, 1L)
  as.integer(sub(pattern, "\\1", hit))
}

test_that("R CVKS bandwidth-method constant matches the C header", {
  ns <- asNamespace("np")
  expect_identical(get("RBWM_CVKS", envir = ns), read_header_define("RBWM_CVKS"))

  body_text <- paste(
    deparse(body(getFromNamespace(".npregbw_call_fixed_degree_core", "np")),
            width.cutoff = 500L),
    collapse = " "
  )
  expect_true(grepl("RBWM_CVKS", body_text, fixed = TRUE))
  expect_false(grepl("identical(objective, \"ks\")) 3L", body_text, fixed = TRUE))
})
