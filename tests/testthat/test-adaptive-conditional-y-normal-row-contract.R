library(npRmpi)

locate_adaptive_y_normal_source <- function() {
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

adaptive_y_normal_body <- function(lines, start_pattern, stop_pattern) {
  start <- grep(start_pattern, lines)
  stop <- grep(stop_pattern, lines)
  expect_length(start, 1L)
  stop <- stop[stop > start][1L]
  expect_true(is.finite(stop))
  paste(lines[start:(stop - 1L)], collapse = "\n")
}

test_that("adaptive ordinary-Gaussian Y normal rows reuse the shared helper", {
  src_file <- locate_adaptive_y_normal_source()
  skip_if(is.null(src_file), "source file src/jksum.c unavailable")
  lines <- readLines(src_file, warn = FALSE)
  body <- adaptive_y_normal_body(
    lines,
    "^static int np_conditional_yrow_from_ctx\\(",
    "^static int np_conditional_y_eval_from_ctx\\("
  )

  expect_match(body, "(BANDWIDTH_den_extern == BW_ADAP_NN)", fixed = TRUE)
  expect_match(body, "(num_var_unordered_extern == 0)", fixed = TRUE)
  expect_match(body, "(num_var_ordered_extern == 0)", fixed = TRUE)
  expect_match(body, "(num_var_continuous_extern > 0)", fixed = TRUE)
  expect_match(body, "(ctx->kernel_cy[0] == 0)", fixed = TRUE)
  expect_match(body, "(ctx->operator_y[0] == OP_NORMAL)", fixed = TRUE)
  expect_match(body, "(!int_cyker_bound_extern)", fixed = TRUE)
  expect_match(body, "(int_TREE_Y != NP_TREE_TRUE)", fixed = TRUE)
  expect_match(body, "np_accel_gauss_adaptive_row_try(", fixed = TRUE)
  expect_match(body, "num_var_continuous_extern,", fixed = TRUE)
  expect_match(body, "num_train,", fixed = TRUE)
  expect_match(body, "1,", fixed = TRUE)
  expect_match(body, "ctx->kw);", fixed = TRUE)
})

test_that("adaptive Y normal admission retains the generic fallback and ownership", {
  src_file <- locate_adaptive_y_normal_source()
  skip_if(is.null(src_file), "source file src/jksum.c unavailable")
  lines <- readLines(src_file, warn = FALSE)
  body <- adaptive_y_normal_body(
    lines,
    "^static int np_conditional_yrow_from_ctx\\(",
    "^static int np_conditional_y_eval_from_ctx\\("
  )

  expect_match(body, "if(!adaptive_gaussian_row &&", fixed = TRUE)
  expect_match(body, "np_conditional_kernel_row(", fixed = TRUE)
  expect_match(body, "np_conditional_push_bounds(", fixed = TRUE)
  expect_match(body, "np_conditional_pop_bounds(", fixed = TRUE)
  expect_match(body, "row_out[orig_j] = ctx->kw[j];", fixed = TRUE)
  expect_false(grepl("malloc|calloc|realloc|alloc_", body))
  expect_false(grepl("OP_CONVOLUTION|OP_INTEGRAL", body))
})
