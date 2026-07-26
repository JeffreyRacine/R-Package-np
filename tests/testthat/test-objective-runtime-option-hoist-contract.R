objective_option_hoist_source <- function() {
  candidates <- c(
    testthat::test_path("..", "..", "src", "jksum.c"),
    testthat::test_path("..", "..", "..", "src", "jksum.c"),
    file.path(Sys.getenv("R_PACKAGE_DIR", ""), "src", "jksum.c"),
    file.path(Sys.getenv("R_PACKAGE_SOURCE", ""), "src", "jksum.c"),
    file.path(getwd(), "src", "jksum.c"),
    file.path(getwd(), "..", "src", "jksum.c")
  )
  candidates <- unique(candidates[nzchar(candidates)])
  hits <- candidates[file.exists(candidates)]
  if (!length(hits))
    return(NULL)
  paste(readLines(hits[[1L]], warn = FALSE), collapse = "\n")
}

test_that("LP objectives freeze the Apple acceleration option only per call", {
  text <- objective_option_hoist_source()
  skip_if(is.null(text), "package C sources unavailable")

  expect_match(text, "int runtime_options_frozen;")
  expect_match(
    text,
    "if\\(\\(outer_pack_ctx == NULL\\) \\|\\|\\s*!outer_pack_ctx->runtime_options_frozen\\)"
  )
  expect_match(text, "const NP_OuterPackCtx frozen_runtime_options")
  expect_match(text, "\\.runtime_options_frozen = 1")
  expect_match(text, "&frozen_runtime_options")
  expect_match(text, "&objective_pack_ctx")
  expect_match(text, "gate_override_ctx,\\s+NULL,\\s+NULL,\\s+0\\);")
})
