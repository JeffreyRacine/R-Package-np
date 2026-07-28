library(npRmpi)

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
  start <- grep(
    "^np_conditional_density_cvls_lp_supertile2_stream\\(",
    lines
  )
  stop <- grep("^int np_conditional_density_cvls_lp_stream\\(", lines)
  expect_length(start, 1L)
  stop <- stop[stop > start][1L]
  expect_true(is.finite(stop))
  paste(lines[start:(stop - 1L)], collapse = "\n")
}

test_that("quadratic output reuses the dead LOO block without changing GEMMs", {
  src_file <- locate_dead_output_source()
  skip_if(is.null(src_file), "source file src/jksum.c unavailable")
  body <- dead_output_body(readLines(src_file, warn = FALSE))

  expect_false(grepl("alloc_vecd\\(block_size\\*block_size\\)", body))
  expect_match(
    body,
    "double * const quad_cross = loo_work[0];",
    fixed = TRUE
  )
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

test_that("dead-output alias begins only after all LOO linear consumers", {
  src_file <- locate_dead_output_source()
  skip_if(is.null(src_file), "source file src/jksum.c unavailable")
  body <- dead_output_body(readLines(src_file, warn = FALSE))

  markers <- c(
    "lin[g] += np_blas_ddot_int",
    "double * const quad_cross = loo_work[0]",
    "np_blas_dgemm_tn_int(ib,",
    "block_terms[block_id[g]]",
    "MPI_Allreduce(&local_fail, &any_fail",
    "MPI_Allreduce(MPI_IN_PLACE"
  )
  positions <- vapply(markers, function(marker) {
    regexpr(marker, body, fixed = TRUE)[[1L]]
  }, integer(1L))
  expect_true(all(positions > 0L))
  expect_true(all(diff(positions) > 0L))
})
