library(np)

locate_dead_output_source <- function() {
  candidates <- c(
    test_path("..", "..", "src", "jksum.c"),
    test_path("..", "..", "..", "src", "jksum.c"),
    file.path(Sys.getenv("R_PACKAGE_DIR", ""), "src", "jksum.c"),
    file.path(Sys.getenv("R_PACKAGE_SOURCE", ""), "src", "jksum.c"),
    file.path(getwd(), "src", "jksum.c"),
    file.path(getwd(), "..", "src", "jksum.c")
  )
  candidates <- unique(candidates[nzchar(candidates)])
  hits <- candidates[file.exists(candidates)]
  if (!length(hits)) NULL else hits[[1L]]
}

dead_output_body <- function(lines) {
  np_test_extract_c_function(
    lines, "np_conditional_density_cvls_lp_supertile2_stream"
  )
}

test_that("quadratic output uses one bounded square tile without extra row slabs", {
  src_file <- locate_dead_output_source()
  skip_if(is.null(src_file), "source file src/jksum.c unavailable")
  body <- dead_output_body(readLines(src_file, warn = FALSE))

  expect_false(grepl("alloc_vecd\\(block_size\\*block_size\\)", body))
  expect_match(
    body,
    "np_cvls_workspace_square_try(block_size, &quad_cross)",
    fixed = TRUE
  )
  expect_false(grepl("loo_work", body, fixed = TRUE))
  expect_false(grepl("full_blocks", body, fixed = TRUE))
  expect_equal(lengths(regmatches(
    body,
    gregexpr("np_blas_dgemm_tn_int\\(", body, perl = TRUE)
  )), 1L)
  expect_match(
    body,
    "np_blas_dgemm_tn_int(ib,",
    fixed = TRUE
  )
})

test_that("bounded square tile preserves LOO consumer order", {
  src_file <- locate_dead_output_source()
  skip_if(is.null(src_file), "source file src/jksum.c unavailable")
  body <- dead_output_body(readLines(src_file, warn = FALSE))

  markers <- c(
    "xblocks[g]) != 0",
    "lin[g] += np_blas_ddot_int",
    "np_blas_dgemm_tn_int(ib,",
    "*cv += quad[g] - 2.0*lin[g]"
  )
  positions <- vapply(markers, function(marker) {
    regexpr(marker, body, fixed = TRUE)[[1L]]
  }, integer(1L))
  expect_true(all(positions > 0L))
  expect_true(all(diff(positions) > 0L))
})
