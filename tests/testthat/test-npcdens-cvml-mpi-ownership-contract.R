npcdens_cvml_mpi_source <- function() {
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

npcdens_cvml_mpi_body <- function(source, start, stop) {
  first <- regexpr(start, source, fixed = TRUE)[[1L]]
  expect_gt(first, 0L)
  remainder <- substr(source, first, nchar(source))
  offset <- regexpr(stop, remainder, fixed = TRUE)[[1L]]
  expect_gt(offset, 0L)
  last <- first + offset - 1L
  substr(source, first, last - 1L)
}

test_that("conditional CVML has one bounded MPI contribution owner", {
  source <- npcdens_cvml_mpi_source()
  skip_if(is.null(source), "package C source unavailable")

  expect_match(source, "np_objective_outer_owned_rows(", fixed = TRUE)
  expect_match(
    source,
    "np_objective_outer_buffer_prepare(",
    fixed = TRUE
  )
  expect_match(
    source,
    "np_objective_outer_buffer_finish(",
    fixed = TRUE
  )
  expect_match(
    source,
    "np_jksum_malloc_array_try(count, sizeof(**buffer))",
    fixed = TRUE
  )
  expect_match(
    source,
    "np_mpi_allreduce_in_place_double(buffer, count, MPI_SUM, label)",
    fixed = TRUE
  )
  expect_false(grepl("alloc_tmatd(num_obs, num_obs)", source, fixed = TRUE))
  expect_false(grepl("diag(num_obs)", source, fixed = TRUE))
})

test_that("prepared outer CVML ownership is kernel-family neutral", {
  source <- npcdens_cvml_mpi_source()
  skip_if(is.null(source), "package C source unavailable")

  policy <- npcdens_cvml_mpi_body(
    source,
    "np_objective_outer_rows_enabled(const int route_ready)",
    "static void np_objective_outer_owned_rows"
  )
  expect_match(policy, "!np_mpi_local_regression_active()", fixed = TRUE)
  expect_match(policy, "route_ready", fixed = TRUE)
  expect_false(grepl("KERNEL_reg_extern", policy, fixed = TRUE))
  expect_false(grepl("KERNEL_den_extern", policy, fixed = TRUE))
  expect_false(grepl("NP_CKERNEL_COORDINATE_CODE", policy, fixed = TRUE))
})

test_that("routed beta CVML uses the canonical prepared outer-row owner", {
  source <- npcdens_cvml_mpi_source()
  skip_if(is.null(source), "package C source unavailable")

  routed_start <- tail(gregexpr(
    "static NP_NOINLINE int np_conditional_density_cvml_continuous_route(",
    source,
    fixed = TRUE
  )[[1L]], 1L)
  routed_stop <- regexpr(
    "typedef struct {\n  int ready;\n  int beta_x;",
    source,
    fixed = TRUE
  )[[1L]]
  expect_gt(routed_start, 0L)
  expect_gt(routed_stop, routed_start)
  routed <- substr(source, routed_start, routed_stop - 1L)
  expect_match(
    routed,
    "np_objective_outer_rows_enabled(int_conditional_prepared_context_extern)",
    fixed = TRUE
  )
  expect_false(grepl("BANDWIDTH_den == BW_FIXED", routed, fixed = TRUE))
  expect_match(
    routed,
    "np_objective_outer_buffer_prepare(",
    fixed = TRUE
  )
  expect_match(routed, "np_objective_outer_owned_rows(", fixed = TRUE)
  expect_match(
    routed,
    "np_objective_outer_buffer_finish(",
    fixed = TRUE
  )
  expect_false(grepl("MPI_Allreduce[[:space:]]*\\(", routed))
  expect_false(grepl("evaluation % iNum_Processors", routed, fixed = TRUE))
})

test_that("all positive-width CVML exits use canonical contribution finish", {
  source <- npcdens_cvml_mpi_source()
  skip_if(is.null(source), "package C source unavailable")

  all_large <- npcdens_cvml_mpi_body(
    source,
    "static int np_conditional_density_cvml_lp_all_large_stream",
    "static int np_conditional_density_cvls_lp_all_large_stream"
  )
  row <- npcdens_cvml_mpi_body(
    source,
    "static int np_conditional_density_cvml_lp_prepared_parallel_stream",
    "int np_conditional_density_cvml_lp_stream"
  )
  block_start <- tail(gregexpr(
    "static int np_conditional_density_cvml_lp_parallel_block_stream",
    source,
    fixed = TRUE
  )[[1L]], 1L)
  block_stop <- regexpr(
    "static int np_conditional_density_cvls_lp_row_stream",
    source,
    fixed = TRUE
  )[[1L]]
  expect_gt(block_start, 0L)
  expect_gt(block_stop, block_start)
  block <- substr(source, block_start, block_stop - 1L)

  for (body in list(all_large, row, block)) {
    expect_match(
      body,
      "np_objective_outer_buffer_finish(",
      fixed = TRUE
    )
    expect_true(
      grepl("use_parallel_rows", body, fixed = TRUE) ||
        grepl("np_objective_outer_owned_rows(i0,\n                                  ib,\n                                  1,", body, fixed = TRUE)
    )
  }

  expect_false(grepl(
    "!int_conditional_prepared_context_extern",
    row,
    fixed = TRUE
  ))
})

test_that("outer CVML ownership suppresses nested kernel parallelism", {
  source <- npcdens_cvml_mpi_source()
  skip_if(is.null(source), "package C source unavailable")

  block_start <- tail(gregexpr(
    "static int np_conditional_density_cvml_lp_parallel_block_stream",
    source,
    fixed = TRUE
  )[[1L]], 1L)
  block_stop <- regexpr(
    "static int np_conditional_density_cvls_lp_row_stream",
    source,
    fixed = TRUE
  )[[1L]]
  block <- substr(source, block_start, block_stop - 1L)

  expect_match(
    block,
    "np_conditional_x_weight_block_stream_core_suppress(",
    fixed = TRUE
  )
  expect_match(
    block,
    "np_conditional_y_block_stream_op_core(vector_scale_factor,",
    fixed = TRUE
  )
  expect_match(block, "OP_NORMAL,\n                                             1,", fixed = TRUE)
  expect_false(grepl(
    "np_conditional_x_weight_block_stream_core(vector_scale_factor,",
    block,
    fixed = TRUE
  ))
})

test_that("noneligible CVML retains the incumbent shared block kernel", {
  source <- npcdens_cvml_mpi_source()
  skip_if(is.null(source), "package C source unavailable")

  shared_start <- tail(gregexpr(
    "static int np_conditional_density_cvml_lp_block_stream(",
    source,
    fixed = TRUE
  )[[1L]], 1L)
  shared_remainder <- substr(source, shared_start, nchar(source))
  shared_length <- regexpr(
    "static int np_conditional_density_cvml_lp_parallel_block_stream(",
    shared_remainder,
    fixed = TRUE
  )[[1L]]
  expect_gt(shared_start, 0L)
  expect_gt(shared_length, 0L)
  shared <- substr(source, shared_start, shared_start + shared_length - 2L)
  expect_match(
    shared,
    "np_conditional_x_weight_block_stream_core(",
    fixed = TRUE
  )
  expect_match(shared, "OP_NORMAL, 0, yblock", fixed = TRUE)
  expect_false(grepl("contributions", shared, fixed = TRUE))
  expect_false(grepl("MPI_Allreduce", shared, fixed = TRUE))
})

test_that("prepared ownership is isolated from the incumbent public route", {
  source <- npcdens_cvml_mpi_source()
  skip_if(is.null(source), "package C source unavailable")

  incumbent <- npcdens_cvml_mpi_body(
    source,
    "int np_conditional_density_cvml_lp_stream(",
    "static int np_conditional_density_cvml_lp_block_stream("
  )
  expect_match(
    incumbent,
    "return np_conditional_density_cvml_lp_prepared_parallel_stream(",
    fixed = TRUE
  )
  expect_match(
    incumbent,
    "!int_conditional_prepared_context_extern;",
    fixed = TRUE
  )
  expect_match(
    incumbent,
    "if(use_parallel_rows && ((i % iNum_Processors) != my_rank))",
    fixed = TRUE
  )
  expect_match(
    incumbent,
    "MPI_Allreduce(&local_cv, cv, 1, MPI_DOUBLE, MPI_SUM, comm[1])",
    fixed = TRUE
  )
  expect_false(grepl("contributions", incumbent, fixed = TRUE))
})

test_that("CVML rank failure uses a common rendezvous", {
  source <- npcdens_cvml_mpi_source()
  skip_if(is.null(source), "package C source unavailable")

  finish <- npcdens_cvml_mpi_body(
    source,
    "np_objective_outer_buffer_finish",
    "static inline int np_regression_cv_scalar_accumulate_scaled_row"
  )
  expect_match(
    source,
    '"NP_RMPI_INJECT_CDEN_CVML_FAIL_RANK"',
    fixed = TRUE
  )
  expect_match(
    finish,
    "MPI_Allreduce(&local_fail, &any_fail, 1, MPI_INT, MPI_MAX, comm[1])",
    fixed = TRUE
  )
  expect_match(finish, "if(any_fail)\n      return 1;", fixed = TRUE)
})
