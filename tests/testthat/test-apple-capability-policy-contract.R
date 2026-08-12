library(np)

locate_apple_policy_source <- function() {
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

apple_policy_source_body <- function(lines, start_pattern, stop_pattern) {
  start <- grep(start_pattern, lines)
  stop <- grep(stop_pattern, lines)
  expect_length(start, 1L)
  expect_gte(length(stop), 1L)
  stop <- stop[stop > start][1L]
  expect_true(is.finite(stop))
  paste(lines[start:(stop - 1L)], collapse = "\n")
}

test_that("Apple-qualified LP selectors own runtime permission", {
  src_file <- locate_apple_policy_source()
  skip_if(is.null(src_file), "source file src/jksum.c unavailable")
  lines <- readLines(src_file, warn = FALSE)
  source <- paste(lines, collapse = "\n")
  dgemv <- apple_policy_source_body(
    lines,
    "^static int np_apple_glp_dgemv_profitable\\(",
    "^static int np_apple_conditional_x_weighted_blas_profitable\\("
  )
  weighted <- apple_policy_source_body(
    lines,
    "^static int np_apple_conditional_x_weighted_blas_profitable\\(",
    "^static int np_reg_alllarge_blas_profitable\\("
  )

  for (body in list(dgemv, weighted)) {
    expect_match(body, "#if NP_ACCEL_GAUSS_COMPILED", fixed = TRUE)
    expect_match(body, "np_mseries_accelerate_enabled_cache", fixed = TRUE)
    expect_false(grepl("np_refresh_mseries_accelerate_option", body,
                       fixed = TRUE))
  }
  expect_false(grepl("np_glp_dgemv_profitable", source, fixed = TRUE))
  expect_false(grepl("np_conditional_x_weighted_blas_profitable", source,
                     fixed = TRUE))
  expect_identical(
    lengths(regmatches(source, gregexpr(
      "np_apple_glp_dgemv_profitable(", source, fixed = TRUE
    ))),
    2L
  )
  expect_identical(
    lengths(regmatches(source, gregexpr(
      "np_apple_conditional_x_weighted_blas_profitable(", source,
      fixed = TRUE
    ))),
    4L
  )
})

test_that("all-large Apple permission is frozen in prepared context", {
  src_file <- locate_apple_policy_source()
  skip_if(is.null(src_file), "source file src/jksum.c unavailable")
  lines <- readLines(src_file, warn = FALSE)
  source <- paste(lines, collapse = "\n")
  row_fit <- apple_policy_source_body(
    lines,
    "^static double np_conditional_lp_all_large_row_fit\\(",
    "^static double np_conditional_lp_all_large_row_fit_basis\\("
  )

  expect_match(source, "int use_apple_dgemv;", fixed = TRUE)
  expect_match(
    source,
    "ctx->use_apple_dgemv =\n    np_apple_glp_dgemv_profitable(",
    fixed = TRUE
  )
  expect_match(row_fit, "if(!ctx->use_apple_dgemv)", fixed = TRUE)
  expect_false(grepl("np_apple_glp_dgemv_profitable", row_fit, fixed = TRUE))
})

test_that("canonical BLAS wrappers remain independent of Apple permission", {
  src_file <- locate_apple_policy_source()
  skip_if(is.null(src_file), "source file src/jksum.c unavailable")
  lines <- readLines(src_file, warn = FALSE)
  body <- apple_policy_source_body(
    lines,
    "^static void np_blas_dgemv_t_int\\(",
    "^static void np_blas_dgemv_n_int\\("
  )

  expect_match(body, "F77_CALL(dgemv)", fixed = TRUE)
  expect_false(grepl("np_mseries_accelerate_enabled_cache", body,
                     fixed = TRUE))
})
