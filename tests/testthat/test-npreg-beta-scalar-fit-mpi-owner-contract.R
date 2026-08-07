npreg_beta_scalar_fit_owner_source <- function() {
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

npreg_beta_scalar_fit_owner_region <- function(source, start, stop) {
  begin <- regexpr(start, source, fixed = TRUE)[[1L]]
  expect_gt(begin, 0L)
  tail <- substr(source, begin, nchar(source))
  finish <- regexpr(stop, tail, fixed = TRUE)[[1L]]
  expect_gt(finish, 0L)
  substr(tail, 1L, finish - 1L)
}

test_that("scalar beta fit keeps one MPI-free mathematical block", {
  source <- npreg_beta_scalar_fit_owner_source()
  skip_if(is.null(source), "package C source unavailable")
  block <- npreg_beta_scalar_fit_owner_region(
    source,
    "static NP_NOINLINE int np_beta_scalar_regression_fit_block_canonical(",
    "#ifdef MPI2\nstatic double **np_beta_scalar_regression_fit_matrix_view("
  )

  expect_match(block, "np_beta_absolute_route(&route_call);", fixed = TRUE)
  expect_match(block, "original_train_is_eval,", fixed = TRUE)
  expect_false(grepl("MPI_", block, fixed = TRUE))
  expect_false(grepl("kernel_estimate_regression_categorical_tree_np(",
                     block, fixed = TRUE))
})

test_that("scalar beta fit owner partitions only complete evaluation rows", {
  source <- npreg_beta_scalar_fit_owner_source()
  skip_if(is.null(source), "package C source unavailable")
  owner <- npreg_beta_scalar_fit_owner_region(
    source,
    "static NP_NOINLINE void np_beta_scalar_regression_fit_canonical(",
    "NP_NOINLINE NP_COLD int np_beta_continuous_bandwidth_prepare_canonical("
  )

  expect_match(
    owner,
    "if(iNum_Processors > 1 && !np_mpi_local_regression_active()) {",
    fixed = TRUE
  )
  expect_match(owner, "np_jksum_mpi_contiguous_row_layout(", fixed = TRUE)
  expect_match(owner, "np_beta_scalar_regression_fit_matrix_view(",
               fixed = TRUE)
  expect_match(owner, "mean + evaluation_start", fixed = TRUE)
  expect_match(owner, "gradient, num_predictors, evaluation_start",
               fixed = TRUE)
  expect_false(grepl("num_obs_train*num_obs_eval", owner, fixed = TRUE))
  expect_false(grepl("num_obs_eval*num_obs_eval", owner, fixed = TRUE))
})

test_that("nearest-neighbor realization is hoisted before rank slicing", {
  source <- npreg_beta_scalar_fit_owner_source()
  skip_if(is.null(source), "package C source unavailable")
  owner <- npreg_beta_scalar_fit_owner_region(
    source,
    "static NP_NOINLINE void np_beta_scalar_regression_fit_canonical(",
    "NP_NOINLINE NP_COLD int np_beta_continuous_bandwidth_prepare_canonical("
  )

  prepare <- regexpr(
    "bandwidth_status = np_beta_bandwidth_prepare_matrix(",
    owner,
    fixed = TRUE
  )[[1L]]
  slice <- regexpr(
    "np_beta_scalar_regression_prepared_view_slice(",
    owner,
    fixed = TRUE
  )[[1L]]
  block <- regexpr(
    "status = np_beta_scalar_regression_fit_block_canonical(",
    owner,
    fixed = TRUE
  )[[1L]]
  expect_gt(prepare, 0L)
  expect_gt(slice, prepare)
  expect_gt(block, slice)
  expect_match(owner, "original_train_is_eval,", fixed = TRUE)
  expect_match(source, "slice->evaluation_offset += evaluation_start;",
               fixed = TRUE)
  expect_match(source, "slice->evaluation_count = evaluation_count;",
               fixed = TRUE)
})

test_that("scalar beta fit rendezvous precedes column-wise transport", {
  source <- npreg_beta_scalar_fit_owner_source()
  skip_if(is.null(source), "package C source unavailable")
  owner <- npreg_beta_scalar_fit_owner_region(
    source,
    "static NP_NOINLINE void np_beta_scalar_regression_fit_canonical(",
    "NP_NOINLINE NP_COLD int np_beta_continuous_bandwidth_prepare_canonical("
  )

  rendezvous <- regexpr(
    "MPI_Allreduce(MPI_IN_PLACE, &failing_rank, 1, MPI_INT, MPI_MIN",
    owner,
    fixed = TRUE
  )[[1L]]
  first_gather <- regexpr(
    "beta scalar regression mean MPI_Allgatherv",
    owner,
    fixed = TRUE
  )[[1L]]
  expect_gt(rendezvous, 0L)
  expect_gt(first_gather, rendezvous)
  for (label in c(
    "beta scalar regression mean MPI_Allgatherv",
    "beta scalar regression standard error MPI_Allgatherv",
    "beta scalar regression gradient MPI_Allgatherv",
    "beta scalar regression gradient standard error MPI_Allgatherv"
  ))
    expect_match(owner, label, fixed = TRUE)
  expect_match(owner, "MPI_Bcast(payload, 7, MPI_INT, failing_rank",
               fixed = TRUE)
})
