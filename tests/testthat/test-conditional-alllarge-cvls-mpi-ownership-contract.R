locate_conditional_alllarge_source <- function() {
  roots <- unique(c(
    Sys.getenv("NP_SOURCE_ROOT", unset = ""),
    normalizePath(file.path(getwd(), "..", ".."), mustWork = FALSE),
    normalizePath(getwd(), mustWork = FALSE)
  ))
  roots <- roots[nzchar(roots)]
  for (root in roots) {
    path <- file.path(root, "src", "jksum.c")
    if (file.exists(path))
      return(path)
  }
  NULL
}

extract_conditional_alllarge_region <- function(source, start, finish) {
  first <- regexpr(start, source, fixed = TRUE)[[1L]]
  tail <- substring(source, first)
  last <- regexpr(finish, tail, fixed = TRUE)[[1L]]

  expect_gt(first, 0L)
  expect_gt(last, 0L)
  substring(tail, 1L, last - 1L)
}

test_that("local density all-large CVLS remains ownership-free", {
  path <- locate_conditional_alllarge_source()
  skip_if(is.null(path), "package source unavailable")
  source <- paste(readLines(path, warn = FALSE), collapse = "\n")
  local <- extract_conditional_alllarge_region(
    source,
    "static int np_conditional_density_cvls_lp_all_large_stream(",
    "static int np_conditional_density_cvls_lp_all_large_parallel_stream("
  )

  expect_match(local, "np_conditional_lp_all_large_build_conv_quad(",
               fixed = TRUE)
  expect_match(local, "*cv += quad - 2.0*lin;", fixed = TRUE)
  expect_false(grepl("np_objective_outer_", local, fixed = TRUE))
  expect_false(grepl("num_train*num_train", local, fixed = TRUE))
})

test_that("MPI density all-large CVLS owns complete rows", {
  path <- locate_conditional_alllarge_source()
  skip_if(is.null(path), "package source unavailable")
  source <- paste(readLines(path, warn = FALSE), collapse = "\n")
  parallel <- extract_conditional_alllarge_region(
    source,
    "static int np_conditional_density_cvls_lp_all_large_parallel_stream(",
    "static int np_conditional_density_cvls_lp_all_large_dispatch("
  )

  expect_match(parallel, "np_objective_outer_preflight_failed(", fixed = TRUE)
  expect_match(parallel, "np_objective_outer_buffer_prepare(", fixed = TRUE)
  expect_match(parallel, "np_objective_outer_owned_rows(", fixed = TRUE)
  expect_match(parallel, "contributions[i] = quad - 2.0*lin;", fixed = TRUE)
  expect_match(parallel, "np_objective_outer_buffer_finish(", fixed = TRUE)
  expect_match(parallel, '"NP_RMPI_INJECT_CDEN_CVLS_FAIL_RANK"',
               fixed = TRUE)
  expect_match(parallel, "np_conditional_y_row_stream_op_core_suppress(",
               fixed = TRUE)
  expect_false(grepl("num_train*num_train", parallel, fixed = TRUE))
  expect_false(grepl("alloc_matd(ctx.num_train, ctx.num_train)",
                     parallel, fixed = TRUE))
})

test_that("density all-large dispatcher isolates local and MPI owners", {
  path <- locate_conditional_alllarge_source()
  skip_if(is.null(path), "package source unavailable")
  source <- paste(readLines(path, warn = FALSE), collapse = "\n")
  dispatch <- extract_conditional_alllarge_region(
    source,
    "static int np_conditional_density_cvls_lp_all_large_dispatch(",
    "static int np_conditional_distribution_cvls_lp_all_large_stream("
  )

  expect_match(dispatch, "np_objective_outer_rows_enabled(", fixed = TRUE)
  expect_match(
    dispatch,
    "np_conditional_density_cvls_lp_all_large_parallel_stream(",
    fixed = TRUE
  )
  expect_match(
    dispatch,
    "np_conditional_density_cvls_lp_all_large_stream(",
    fixed = TRUE
  )
  expect_match(
    source,
    "np_conditional_density_cvls_lp_all_large_dispatch(vector_scale_factor, cv)",
    fixed = TRUE
  )
})
