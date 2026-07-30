library(np)

locate_cvls_workspace_source <- function(filename) {
  candidates <- c(
    test_path("..", "..", "src", filename),
    test_path("..", "..", "..", "src", filename),
    file.path(Sys.getenv("R_PACKAGE_DIR", ""), "src", filename),
    file.path(Sys.getenv("R_PACKAGE_SOURCE", ""), "src", filename),
    file.path(getwd(), "src", filename),
    file.path(getwd(), "..", "src", filename)
  )
  candidates <- unique(candidates[nzchar(candidates)])
  hits <- candidates[file.exists(candidates)]
  if (!length(hits)) NULL else hits[[1L]]
}

test_that("conditional LP tile planner owns a truthful total budget", {
  src_file <- locate_cvls_workspace_source("jksum.c")
  skip_if(is.null(src_file), "source file src/jksum.c unavailable")
  src <- paste(readLines(src_file, warn = FALSE), collapse = "\n")

  expect_match(
    src,
    "#define NP_CONDITIONAL_LP_TILE_BUDGET_BYTES",
    fixed = TRUE
  )
  expect_match(
    src,
    "((size_t)512 * (size_t)1024 * (size_t)1024)",
    fixed = TRUE
  )
  expect_match(src, "np_native_bounded_tile_width(", fixed = TRUE)
  expect_match(
    src,
    "np_conditional_lp_cvls_preferred_block_size(n)",
    fixed = TRUE
  )
  expect_match(
    src,
    "np_conditional_lp_cvls_block_size(num_obs, 6U, 0U)",
    fixed = TRUE
  )
  expect_match(
    src,
    "np_conditional_lp_cvls_block_size(num_train, 5U, 1U)",
    fixed = TRUE
  )
  expect_match(
    src,
    "np_conditional_lp_cvls_block_size(num_obs, 2U, 0U)",
    fixed = TRUE
  )
})

test_that("all conditional LP acceleration slabs use typed allocation", {
  src_file <- locate_cvls_workspace_source("jksum.c")
  skip_if(is.null(src_file), "source file src/jksum.c unavailable")
  lines <- readLines(src_file, warn = FALSE)
  start <- grep(
    "^static int np_conditional_density_cvml_lp_block_stream\\(",
    lines
  )
  stop <- grep("^#undef NP_CDENS_ADAP_WIDTH4_NOINLINE$", lines)
  expect_gte(length(start), 1L)
  start <- tail(start, 1L)
  stop <- stop[stop > start][1L]
  expect_true(is.finite(stop))
  body <- paste(lines[start:stop], collapse = "\n")

  expect_match(body, "np_cvls_workspace_matrix_try(", fixed = TRUE)
  expect_match(body, "np_cvls_workspace_vector_try(", fixed = TRUE)
  expect_match(body, "np_cvls_workspace_square_try(", fixed = TRUE)
  expect_false(grepl(
    "(size_t)block_size*(size_t)block_size",
    body,
    fixed = TRUE
  ))
  expect_false(grepl("np_optional_tmatd", body, fixed = TRUE))
  expect_false(grepl(
    "alloc_tmatd\\(num_(obs|train), block_size\\)",
    body,
    perl = TRUE
  ))
  expect_false(grepl(
    "alloc_vecd\\(block_size\\*block_size\\)",
    body,
    perl = TRUE
  ))
  expect_match(
    body,
    "workspace_status == NP_CVLS_WORKSPACE_UNAVAILABLE",
    fixed = TRUE
  )
  expect_match(
    body,
    "return np_conditional_density_cvls_lp_row_stream",
    fixed = TRUE
  )
  expect_match(
    body,
    "return np_conditional_distribution_cvls_lp_row_stream",
    fixed = TRUE
  )
  expect_match(
    body,
    "np_size_mul_checked((size_t)width, (size_t)width, &count)",
    fixed = TRUE
  )
})

test_that("native tile accounting includes data, pointers, and BLAS output", {
  header_file <- locate_cvls_workspace_source("np_native_safety.h")
  skip_if(is.null(header_file), "source file src/np_native_safety.h unavailable")
  header <- paste(readLines(header_file, warn = FALSE), collapse = "\n")

  expect_match(header, "np_native_tile_workspace_bytes(", fixed = TRUE)
  expect_match(header, "linear_matrix_count", fixed = TRUE)
  expect_match(header, "square_matrix_count", fixed = TRUE)
  expect_match(header, "sizeof(double *)", fixed = TRUE)
  expect_match(header, "np_native_bounded_tile_width(", fixed = TRUE)
  expect_match(header, "for(width = preferred_width; width > 0; width /= 2)",
               fixed = TRUE)
  expect_match(header, "np_native_malloc_column_matrix(", fixed = TRUE)
  expect_match(header, "NP_NATIVE_ALLOC_UNAVAILABLE", fixed = TRUE)
})
