beta_npksum_owner_source <- function() {
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

beta_npksum_owner_region <- function(source, start, stop) {
  begin <- regexpr(start, source, fixed = TRUE)[[1L]]
  expect_gt(begin, 0L)
  tail <- substr(source, begin, nchar(source))
  finish <- regexpr(stop, tail, fixed = TRUE)[[1L]]
  expect_gt(finish, 0L)
  substr(tail, 1L, finish - 1L)
}

test_that("absolute beta mathematics remains MPI-free", {
  source <- beta_npksum_owner_source()
  skip_if(is.null(source), "package C source unavailable")
  local_route <- beta_npksum_owner_region(
    source,
    "static int np_beta_absolute_route_body(",
    "#ifdef MPI2\nstatic double **np_beta_absolute_route_matrix_view("
  )

  expect_match(local_route, "np_continuous_kernel_beta_factor_row(",
               fixed = TRUE)
  expect_false(grepl("MPI_", local_route, fixed = TRUE))
  expect_false(grepl("np_mpi_", local_route, fixed = TRUE))
})

test_that("public beta owner admits fixed and GNN but not ANN", {
  source <- beta_npksum_owner_source()
  skip_if(is.null(source), "package C source unavailable")
  dispatcher <- beta_npksum_owner_region(
    source,
    "static int np_beta_absolute_route_dispatch(",
    "/*\n  Keep this adjacent hot engine"
  )

  expect_match(dispatcher, "!suppress_parallel", fixed = TRUE)
  expect_match(dispatcher, "!np_mpi_local_regression_active()", fixed = TRUE)
  expect_match(dispatcher, "call->bandwidth_mode == BW_FIXED", fixed = TRUE)
  expect_match(dispatcher, "call->bandwidth_mode == BW_GEN_NN", fixed = TRUE)
  expect_false(grepl("call->bandwidth_mode == BW_ADAP_NN",
                     dispatcher, fixed = TRUE))
})

test_that("beta owner transports only complete row-major outputs", {
  source <- beta_npksum_owner_source()
  skip_if(is.null(source), "package C source unavailable")
  gather <- beta_npksum_owner_region(
    source,
    "static int np_beta_absolute_route_gather_rows(",
    "/*\n * Coarse operation owner"
  )
  owner <- beta_npksum_owner_region(
    source,
    "static int np_beta_absolute_route_mpi_owner(",
    "#endif\n\nstatic int np_beta_absolute_route_dispatch("
  )

  expect_match(gather, "np_jksum_mpi_contiguous_row_layout(", fixed = TRUE)
  expect_match(owner, "np_beta_absolute_route_gather_rows(", fixed = TRUE)
  expect_match(owner, "local_call.num_obs_eval = evaluation_count;",
               fixed = TRUE)
  expect_match(owner, "local_call.leave_one_out_offset += evaluation_start;",
               fixed = TRUE)
  expect_match(owner, "local_call.leave_one_out_train_is_eval =",
               fixed = TRUE)
  expect_match(owner, '"NP_RMPI_INJECT_BETA_NPKSUM_FAIL_RANK"',
               fixed = TRUE)
  for (label in c(
    "beta kernel sum MPI_Allgatherv",
    "beta dual-power kernel sum MPI_Allgatherv",
    "beta centered moment MPI_Allgatherv",
    "beta kernel weights MPI_Allgatherv"
  ))
    expect_match(owner, label, fixed = TRUE)
  expect_false(grepl("num_obs_eval * num_obs_eval", owner, fixed = TRUE))
  expect_false(grepl("num_obs_train * num_obs_eval", owner, fixed = TRUE))
})

test_that("beta owner rendezvous precedes output transport", {
  source <- beta_npksum_owner_source()
  skip_if(is.null(source), "package C source unavailable")
  owner <- beta_npksum_owner_region(
    source,
    "static int np_beta_absolute_route_mpi_owner(",
    "#endif\n\nstatic int np_beta_absolute_route_dispatch("
  )
  rendezvous <- regexpr(
    "MPI_Allreduce(MPI_IN_PLACE, &failing_rank", owner, fixed = TRUE
  )[[1L]]
  first_gather <- regexpr(
    "beta kernel sum MPI_Allgatherv", owner, fixed = TRUE
  )[[1L]]
  diagnostic_reduce <- regexpr(
    "MPI_Allreduce(MPI_IN_PLACE, &undefined_count", owner, fixed = TRUE
  )[[1L]]

  expect_gt(rendezvous, 0L)
  expect_gt(first_gather, rendezvous)
  expect_gt(diagnostic_reduce, first_gather)
  expect_match(owner, "MPI_Bcast(failure_payload, 5, MPI_INT",
               fixed = TRUE)
})
