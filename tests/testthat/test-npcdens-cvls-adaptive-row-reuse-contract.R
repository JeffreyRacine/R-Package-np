library(npRmpi)

locate_adaptive_row_reuse_sources <- function() {
  roots <- c(
    test_path("..", ".."),
    test_path("..", "..", ".."),
    Sys.getenv("R_PACKAGE_DIR", ""),
    Sys.getenv("R_PACKAGE_SOURCE", ""),
    getwd(),
    file.path(getwd(), "..")
  )
  roots <- unique(roots[nzchar(roots)])
  for (root in roots) {
    files <- file.path(root, "src", c(
      "jksum.c", "jksum_lp_row.c", "jksum_lp_row.h"
    ))
    if (all(file.exists(files)))
      return(files)
  }
  NULL
}

adaptive_source_body <- function(lines, start_pattern, stop_pattern) {
  start <- grep(start_pattern, lines)
  stop <- grep(stop_pattern, lines)
  expect_length(start, 1L)
  expect_gte(length(stop), 1L)
  stop <- stop[stop > start][1L]
  expect_true(is.finite(stop))
  paste(lines[start:(stop - 1L)], collapse = "\n")
}

test_that("adaptive CVLS consumes only canonical delete-one context rows", {
  files <- locate_adaptive_row_reuse_sources()
  skip_if(is.null(files), "package source files unavailable")
  jksum <- readLines(files[[1L]], warn = FALSE)
  source <- paste(jksum, collapse = "\n")

  expect_false(grepl("np_lp_delete_smoother_row", source, fixed = TRUE))
  expect_false(grepl("np_conditional_xrow_full_from_ctx", source, fixed = TRUE))

  row_body <- adaptive_source_body(
    jksum,
    "^static int np_conditional_density_cvls_lp_row_stream\\(",
    "^static int np_conditional_density_cvls_lp_adap_block_stream\\("
  )
  block_body <- adaptive_source_body(
    jksum,
    "^static int np_conditional_density_cvls_lp_adap_block_stream\\(",
    "^np_conditional_density_cvls_categorical_profile_stream\\("
  )

  for (body in list(row_body, block_body)) {
    expect_match(body, "np_conditional_xrow_from_ctx", fixed = TRUE)
    expect_false(grepl("np_conditional_xrow_from_ctx_impl", body, fixed = TRUE))
    expect_false(grepl("np_lp_delete_smoother_row", body, fixed = TRUE))
  }
})

test_that("adaptive row context owns signed delete-one conversion", {
  files <- locate_adaptive_row_reuse_sources()
  skip_if(is.null(files), "package source files unavailable")
  lines <- readLines(files[[1L]], warn = FALSE)
  implementation <- npRmpi_test_extract_c_function(
    lines, "np_conditional_xrow_influence"
  )
  row_source <- paste(readLines(files[[2L]], warn = FALSE), collapse = "\n")
  header <- paste(readLines(files[[3L]], warn = FALSE), collapse = "\n")

  expect_match(implementation, "if(drop_eval_self){", fixed = TRUE)
  expect_match(
    implementation,
    "if(!np_lp_delete_denominator(row_out[eval_idx], &denominator))",
    fixed = TRUE
  )
  expect_match(
    implementation,
    "j == eval_pos ? 0.0 : row_out[orig_j]/denominator",
    fixed = TRUE
  )
  expect_false(grepl("NZD", implementation, fixed = TRUE))
  expect_false(grepl("np_lp_delete_smoother_row", row_source, fixed = TRUE))
  expect_false(grepl("np_lp_delete_smoother_row", header, fixed = TRUE))
})
