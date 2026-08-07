locate_conditional_profile_source_root <- function() {
  roots <- c(
    testthat::test_path("..", ".."),
    testthat::test_path("..", "..", ".."),
    Sys.getenv("R_PACKAGE_SOURCE", unset = ""),
    getwd()
  )
  roots <- roots[
    file.exists(file.path(roots, "src", "jksum.c")) &
      file.exists(file.path(roots, "src", "np.c"))
  ]
  if (!length(roots)) NULL else roots[[1L]]
}

extract_profile_implementation <- function(source, start_pattern, end_pattern) {
  starts <- gregexpr(start_pattern, source, fixed = TRUE)[[1L]]
  start <- tail(starts[starts > 0L], 1L)
  if (!length(start)) return(NULL)
  remainder <- substring(source, start)
  finish <- regexpr(end_pattern, remainder, fixed = TRUE)[[1L]]
  if (finish <= 0L) return(NULL)
  substring(remainder, 1L, finish - 1L)
}

test_that("conditional profile indexes belong to prepared objective state", {
  root <- locate_conditional_profile_source_root()
  skip_if(is.null(root), "package source unavailable")

  jksum <- paste(
    readLines(file.path(root, "src", "jksum.c"), warn = FALSE),
    collapse = "\n"
  )
  np <- paste(
    readLines(file.path(root, "src", "np.c"), warn = FALSE),
    collapse = "\n"
  )
  headers <- paste(
    readLines(file.path(root, "src", "headers.h"), warn = FALSE),
    collapse = "\n"
  )

  density <- extract_profile_implementation(
    jksum,
    paste0(
      "static NPConditionalProfileCvStatus\n",
      "np_conditional_density_cvls_categorical_profile_stream("
    ),
    paste0(
      "\nstatic NPConditionalProfileCvStatus\n",
      "np_conditional_distribution_cvls_categorical_profile_stream("
    )
  )
  distribution <- extract_profile_implementation(
    jksum,
    paste0(
      "static NPConditionalProfileCvStatus\n",
      "np_conditional_distribution_cvls_categorical_profile_stream("
    ),
    "\nint np_conditional_distribution_cvls_lp_stream("
  )
  expect_false(is.null(density))
  expect_false(is.null(distribution))

  expect_match(jksum, "NPConditionalProfileIndexState", fixed = TRUE)
  expect_match(
    jksum,
    "np_conditional_profile_index_state_matches",
    fixed = TRUE
  )
  expect_match(
    jksum,
    "np_conditional_profile_index_cache_clear_extern",
    fixed = TRUE
  )
  expect_match(headers,
               "np_conditional_profile_index_cache_clear_extern",
               fixed = TRUE)

  for (implementation in list(density, distribution)) {
    expect_match(
      implementation,
      "np_conditional_profile_index_state_prepare",
      fixed = TRUE
    )
    expect_match(
      implementation,
      "int_conditional_prepared_context_extern",
      fixed = TRUE
    )
    expect_false(grepl(
      "np_build_discrete_profile_index",
      implementation,
      fixed = TRUE
    ))
    expect_false(grepl(
      "alloc_vecd(nprof_xy*nprof_xy)",
      implementation,
      fixed = TRUE
    ))
  }

  density_starts <- gregexpr(
    "static void np_conditional_density_prepared_context_clear_internal(void)",
    np,
    fixed = TRUE
  )[[1L]]
  density_start <- tail(density_starts[density_starts > 0L], 1L)
  density_tail <- substring(np, density_start)
  density_clear <- regexpr(
    "np_conditional_profile_index_cache_clear_extern();",
    density_tail,
    fixed = TRUE
  )[[1L]]
  density_free <- regexpr(
    "free_mat(matrix_X_unordered_train_extern",
    density_tail,
    fixed = TRUE
  )[[1L]]
  expect_gt(density_start, 0L)
  expect_gt(density_clear, 0L)
  expect_gt(density_free, density_clear)

  distribution_starts <- gregexpr(
    "static void np_conditional_distribution_prepared_context_destroy(",
    np,
    fixed = TRUE
  )[[1L]]
  distribution_start <- tail(
    distribution_starts[distribution_starts > 0L],
    1L
  )
  distribution_tail <- substring(np, distribution_start)
  distribution_clear <- regexpr(
    "np_conditional_profile_index_cache_clear_extern();",
    distribution_tail,
    fixed = TRUE
  )[[1L]]
  distribution_free <- regexpr(
    "free_mat(context->matrix_y_unordered_train",
    distribution_tail,
    fixed = TRUE
  )[[1L]]
  expect_gt(distribution_start, 0L)
  expect_gt(distribution_clear, 0L)
  expect_gt(distribution_free, distribution_clear)
  expect_match(
    distribution_tail,
    "int_conditional_prepared_context_extern = 0;",
    fixed = TRUE
  )
  expect_match(
    np,
    paste0(
      "prepared_context->active = 1;\n",
      "  int_conditional_prepared_context_extern = 1;"
    ),
    fixed = TRUE
  )
})
