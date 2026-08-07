npudens_cvml_mpi_source <- function() {
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

npudens_cvml_mpi_body <- function(source, start, stop) {
  first <- regexpr(start, source, fixed = TRUE)[[1L]]
  expect_gt(first, 0L)
  remainder <- substr(source, first, nchar(source))
  offset <- regexpr(stop, remainder, fixed = TRUE)[[1L]]
  expect_gt(offset, 0L)
  substr(source, first, first + offset - 2L)
}

test_that("unconditional beta CVML owns complete evaluation rows", {
  source <- npudens_cvml_mpi_source()
  skip_if(is.null(source), "package C source unavailable")

  body <- npudens_cvml_mpi_body(
    source,
    "static int np_density_cvml_beta_route(",
    "static int np_density_cvls_beta_cross_term("
  )
  expect_match(body, "np_objective_outer_rows_enabled(1)", fixed = TRUE)
  expect_match(body, "np_objective_outer_owned_rows(", fixed = TRUE)
  expect_match(body, "np_objective_outer_buffer_prepare(", fixed = TRUE)
  expect_match(body, "np_objective_outer_buffer_finish(", fixed = TRUE)
  expect_match(
    body,
    '"NP_RMPI_INJECT_UDEN_CVML_FAIL_RANK"',
    fixed = TRUE
  )
  expect_false(grepl("evaluation % iNum_Processors", body, fixed = TRUE))
  expect_false(grepl("MPI_Allreduce(MPI_IN_PLACE, cv", body, fixed = TRUE))
})

test_that("unconditional beta CVML has one shared row contribution finisher", {
  source <- npudens_cvml_mpi_source()
  skip_if(is.null(source), "package C source unavailable")

  finisher <- npudens_cvml_mpi_body(
    source,
    "static inline int np_density_cvml_beta_row_contribution(",
    "static int np_density_cvml_beta_route("
  )
  body <- npudens_cvml_mpi_body(
    source,
    "static int np_density_cvml_beta_route(",
    "static int np_density_cvls_beta_cross_term("
  )
  expect_match(
    finisher,
    "np_beta_scaled_row_context_fill_omitting(",
    fixed = TRUE
  )
  expect_match(
    finisher,
    "np_guarded_cvml_log_contribution(",
    fixed = TRUE
  )
  expect_equal(
    lengths(regmatches(body, gregexpr(
      "np_density_cvml_beta_row_contribution(",
      body,
      fixed = TRUE
    ))),
    1L
  )
})

test_that("unconditional beta CVML failure rendezvous precedes result transport", {
  source <- npudens_cvml_mpi_source()
  skip_if(is.null(source), "package C source unavailable")

  body <- npudens_cvml_mpi_body(
    source,
    "static int np_density_cvml_beta_route(",
    "static int np_density_cvls_beta_cross_term("
  )
  prepare_fail <- regexpr(
    "MPI_Allreduce(&local_fail, &any_fail, 1, MPI_INT, MPI_MAX, comm[1])",
    body,
    fixed = TRUE
  )[[1L]]
  finish <- regexpr("np_objective_outer_buffer_finish(", body, fixed = TRUE)[[1L]]
  expect_gt(prepare_fail, 0L)
  expect_gt(finish, prepare_fail)
  expect_match(body, "if(any_fail)\n      goto cleanup;", fixed = TRUE)
})

test_that("unconditional beta CVML adds only linear contribution storage", {
  source <- npudens_cvml_mpi_source()
  skip_if(is.null(source), "package C source unavailable")

  body <- npudens_cvml_mpi_body(
    source,
    "static int np_density_cvml_beta_route(",
    "static int np_density_cvls_beta_cross_term("
  )
  expect_match(body, "(size_t)num_obs", fixed = TRUE)
  expect_false(grepl("num_obs * num_obs", body, fixed = TRUE))
  expect_false(grepl("num_obs\\*num_obs", body))
})
