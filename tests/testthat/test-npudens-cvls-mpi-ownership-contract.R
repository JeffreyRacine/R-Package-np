npudens_cvls_mpi_source <- function() {
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

npudens_cvls_mpi_body <- function(source, start, stop) {
  first <- regexpr(start, source, fixed = TRUE)[[1L]]
  expect_gt(first, 0L)
  remainder <- substr(source, first, nchar(source))
  offset <- regexpr(stop, remainder, fixed = TRUE)[[1L]]
  expect_gt(offset, 0L)
  substr(source, first, first + offset - 2L)
}

test_that("unconditional beta CVLS cross term owns complete delete-one rows", {
  source <- npudens_cvls_mpi_source()
  skip_if(is.null(source), "package C source unavailable")

  local_body <- npudens_cvls_mpi_body(
    source,
    "static int np_density_cvls_beta_cross_term(",
    "static int np_density_cvls_beta_cross_term_mpi("
  )
  mpi_body <- npudens_cvls_mpi_body(
    source,
    "static int np_density_cvls_beta_cross_term_mpi(",
    "int np_kernel_estimate_density_categorical_leave_one_out_cv("
  )
  expect_false(grepl("np_objective_outer_", local_body, fixed = TRUE))
  expect_false(grepl("MPI_", local_body, fixed = TRUE))
  expect_match(mpi_body, "np_objective_outer_rows_enabled(1)", fixed = TRUE)
  expect_match(mpi_body, "np_objective_outer_owned_rows(", fixed = TRUE)
  expect_match(mpi_body, "np_objective_outer_buffer_prepare(", fixed = TRUE)
  expect_match(mpi_body, "np_objective_outer_buffer_finish(", fixed = TRUE)
  expect_match(
    mpi_body,
    '"NP_RMPI_INJECT_UDEN_CVLS_CROSS_FAIL_RANK"',
    fixed = TRUE
  )
  expect_match(mpi_body, "contributions[evaluation]", fixed = TRUE)
  expect_false(grepl("evaluation % iNum_Processors", mpi_body, fixed = TRUE))
  expect_false(grepl("MPI_Allreduce(MPI_IN_PLACE, cross_term", mpi_body,
                     fixed = TRUE))
})

test_that("unconditional beta CVLS cross storage remains linear", {
  source <- npudens_cvls_mpi_source()
  skip_if(is.null(source), "package C source unavailable")

  body <- npudens_cvls_mpi_body(
    source,
    "static int np_density_cvls_beta_cross_term_mpi(",
    "int np_kernel_estimate_density_categorical_leave_one_out_cv("
  )
  expect_match(body, "(size_t)num_obs", fixed = TRUE)
  expect_false(grepl("num_obs * num_obs", body, fixed = TRUE))
  expect_false(grepl("num_obs\\*num_obs", body))
})

test_that("unconditional beta CVLS selects one terminal cross-term sibling", {
  source <- npudens_cvls_mpi_source()
  skip_if(is.null(source), "package C source unavailable")

  body <- npudens_cvls_mpi_body(
    source,
    "int np_kernel_estimate_density_categorical_convolution_cv(",
    "void kernel_estimate_dens_dist_categorical_np("
  )
  expect_match(
    body,
    "if(np_objective_outer_rows_enabled(1)) {",
    fixed = TRUE
  )
  expect_equal(
    lengths(regmatches(body, gregexpr(
      "np_density_cvls_beta_cross_term_mpi(",
      body,
      fixed = TRUE
    ))),
    1L
  )
  expect_equal(
    lengths(regmatches(body, gregexpr(
      "np_density_cvls_beta_cross_term(",
      body,
      fixed = TRUE
    ))),
    1L
  )
})

test_that("unconditional beta CVLS quadrature owns bounded row tiles", {
  source <- npudens_cvls_mpi_source()
  skip_if(is.null(source), "package C source unavailable")

  body <- npudens_cvls_mpi_body(
    source,
    "static int np_density_cvls_bounded_i1_quadrature_general_beta_mpi(",
    "int np_kernel_estimate_density_categorical_convolution_cv("
  )
  expect_match(body, "np_objective_outer_rows_enabled(1)", fixed = TRUE)
  expect_match(body, "np_objective_outer_owned_rows(", fixed = TRUE)
  expect_match(body, "np_objective_outer_buffer_prepare(", fixed = TRUE)
  expect_match(body, "(size_t)block_size", fixed = TRUE)
  expect_match(body, "np_objective_outer_buffer_finish(", fixed = TRUE)
  expect_match(
    body,
    '"NP_RMPI_INJECT_UDEN_CVLS_QUAD_FAIL_RANK"',
    fixed = TRUE
  )
  expect_match(
    body,
    "eval_start + (size_t)owned_offset",
    fixed = TRUE
  )
  expect_match(
    body,
    "nord,\n                               1,\n                               vector_scale_factor",
    fixed = TRUE
  )
  expect_false(grepl("(size_t)total_eval,\n", body, fixed = TRUE))
  expect_false(grepl("total_eval * total_eval", body, fixed = TRUE))
  expect_false(grepl("num_obs * num_obs", body, fixed = TRUE))
})

test_that("unconditional beta CVLS selects one terminal quadrature sibling", {
  source <- npudens_cvls_mpi_source()
  skip_if(is.null(source), "package C source unavailable")

  body <- npudens_cvls_mpi_body(
    source,
    "int np_kernel_estimate_density_categorical_convolution_cv(",
    "void kernel_estimate_dens_dist_categorical_np("
  )
  expect_match(
    body,
    "exact_beta_route && np_objective_outer_rows_enabled(1)",
    fixed = TRUE
  )
  expect_equal(
    lengths(regmatches(body, gregexpr(
      "np_density_cvls_bounded_i1_quadrature_general_beta_mpi(",
      body,
      fixed = TRUE
    ))),
    1L
  )
  expect_equal(
    lengths(regmatches(body, gregexpr(
      "np_density_cvls_bounded_i1_quadrature_general(",
      body,
      fixed = TRUE
    ))),
    1L
  )
})
