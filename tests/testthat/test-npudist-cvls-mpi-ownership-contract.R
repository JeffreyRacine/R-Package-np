npudist_cvls_mpi_source <- function() {
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
  if (!length(hits)) return(NULL)
  paste(readLines(hits[[1L]], warn = FALSE), collapse = "\n")
}

npudist_cvls_mpi_body <- function(source, start, stop) {
  first <- regexpr(start, source, fixed = TRUE)[[1L]]
  expect_gt(first, 0L)
  remainder <- substr(source, first, nchar(source))
  offset <- regexpr(stop, remainder, fixed = TRUE)[[1L]]
  expect_gt(offset, 0L)
  substr(source, first, first + offset - 2L)
}

test_that("unconditional beta distribution CVLS owns complete evaluations", {
  source <- npudist_cvls_mpi_source()
  skip_if(is.null(source), "package C source unavailable")

  body <- npudist_cvls_mpi_body(
    source,
    "static int np_distribution_cvls_continuous_route_mpi(",
    "double np_kernel_estimate_distribution_ls_cv("
  )
  expect_match(body, "np_objective_outer_rows_enabled(1)", fixed = TRUE)
  expect_match(body, "np_objective_outer_owned_rows(", fixed = TRUE)
  expect_match(body, "np_objective_outer_buffer_prepare(", fixed = TRUE)
  expect_match(body, "np_objective_outer_buffer_finish(", fixed = TRUE)
  expect_match(
    body,
    '"NP_RMPI_INJECT_UDIST_CVLS_FAIL_RANK"',
    fixed = TRUE
  )
  expect_match(
    body,
    "evaluation, 0, num_obs_train - 1, cdfontrain",
    fixed = TRUE
  )
  expect_match(body, "contributions[evaluation]", fixed = TRUE)
  expect_false(grepl("evaluation % iNum_Processors", body, fixed = TRUE))
  expect_false(grepl("MPI_Allreduce(MPI_IN_PLACE, cv", body,
                     fixed = TRUE))
})

test_that("unconditional beta distribution CVLS storage is linear", {
  source <- npudist_cvls_mpi_source()
  skip_if(is.null(source), "package C source unavailable")

  body <- npudist_cvls_mpi_body(
    source,
    "static int np_distribution_cvls_continuous_route_mpi(",
    "double np_kernel_estimate_distribution_ls_cv("
  )
  expect_match(body, "(size_t)num_obs_eval", fixed = TRUE)
  expect_match(body, "(size_t)num_obs_train", fixed = TRUE)
  expect_false(grepl("num_obs_train * num_obs_eval", body, fixed = TRUE))
  expect_false(grepl("num_obs_train\\*num_obs_eval", body))
  expect_false(grepl("num_obs_eval * num_obs_eval", body, fixed = TRUE))
})

test_that("unconditional beta distribution CVLS selects one terminal sibling", {
  source <- npudist_cvls_mpi_source()
  skip_if(is.null(source), "package C source unavailable")

  body <- npudist_cvls_mpi_body(
    source,
    "double np_kernel_estimate_distribution_ls_cv(",
    "static int np_conditional_profile_cache_priority_compare("
  )
  expect_match(
    body,
    "if(np_objective_outer_rows_enabled(1)) {",
    fixed = TRUE
  )
  expect_equal(
    lengths(regmatches(body, gregexpr(
      "np_distribution_cvls_continuous_route_mpi(",
      body,
      fixed = TRUE
    ))),
    1L
  )
  expect_equal(
    lengths(regmatches(body, gregexpr(
      "np_distribution_cvls_continuous_route(",
      body,
      fixed = TRUE
    ))),
    1L
  )
})
