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

test_that("adaptive CVLS derives deleted rows only from completed full rows", {
  files <- locate_adaptive_row_reuse_sources()
  skip_if(is.null(files), "package source files unavailable")
  jksum <- readLines(files[[1L]], warn = FALSE)
  source <- paste(jksum, collapse = "\n")

  expect_equal(lengths(regmatches(
    source,
    gregexpr("np_lp_delete_smoother_row\\(", source, perl = TRUE)
  )), 4L)

  row_body <- adaptive_source_body(
    jksum,
    "^static int np_conditional_density_cvls_lp_row_stream\\(",
    "^static int np_conditional_density_cvls_lp_adap_block_stream\\("
  )
  block_body <- adaptive_source_body(
    jksum,
    "^static int np_conditional_density_cvls_lp_adap_block_stream\\(",
    "^static int np_conditional_density_cvls_categorical_profile_stream\\("
  )

  for (body in list(row_body, block_body)) {
    expect_match(
      body,
      "np_lp_engine_extern == NP_LP_ENGINE_GENERAL",
      fixed = TRUE
    )
    full_pos <- regexpr(
      "np_conditional_xrow_full_from_ctx",
      body,
      fixed = TRUE
    )[[1L]]
    delete_pos <- regexpr(
      "np_lp_delete_smoother_row",
      body,
      fixed = TRUE
    )[[1L]]
    expect_true(full_pos > 0L && delete_pos > full_pos)
  }
})

test_that("shared deletion utility is hidden, signed, bounded, and allocation-free", {
  files <- locate_adaptive_row_reuse_sources()
  skip_if(is.null(files), "package source files unavailable")
  implementation_lines <- readLines(files[[2L]], warn = FALSE)
  implementation_start <- grep(
    "^int np_lp_delete_smoother_row\\(",
    implementation_lines
  )
  expect_length(implementation_start, 1L)
  implementation <- paste(
    implementation_lines[implementation_start:length(implementation_lines)],
    collapse = "\n"
  )
  header <- paste(readLines(files[[3L]], warn = FALSE), collapse = "\n")

  expect_match(
    header,
    "attribute_hidden int np_lp_delete_smoother_row",
    fixed = TRUE
  )
  expect_match(
    implementation,
    "if(!np_lp_delete_denominator(full_row[eval_idx], &den))",
    fixed = TRUE
  )
  expect_false(grepl("NZD", implementation, fixed = TRUE))
  expect_match(
    implementation,
    "(j == eval_idx) ? 0.0 : full_row[j]/den",
    fixed = TRUE
  )
  expect_match(implementation, "n <= 0", fixed = TRUE)
  expect_match(implementation, "eval_idx >= n", fixed = TRUE)
  expect_false(grepl("malloc", implementation, fixed = TRUE))
  expect_false(grepl("calloc", implementation, fixed = TRUE))
})
