library(np)

locate_block_plan_source <- function(name) {
  candidates <- c(
    test_path("..", "..", "src", name),
    test_path("..", "..", "..", "src", name),
    file.path(Sys.getenv("R_PACKAGE_DIR", ""), "src", name),
    file.path(Sys.getenv("R_PACKAGE_SOURCE", ""), "src", name),
    file.path(getwd(), "src", name),
    file.path(getwd(), "..", "src", name)
  )
  candidates <- unique(candidates[nzchar(candidates)])
  hits <- candidates[file.exists(candidates)]
  if (length(hits) == 0L) {
    return(NULL)
  }
  hits[[1L]]
}

extract_native_function <- function(lines, start_pattern, stop_pattern) {
  start <- grep(start_pattern, lines)
  stop <- grep(stop_pattern, lines)
  expect_length(start, 1L)
  expect_length(stop, 1L)
  expect_gt(stop, start)
  lines[start:(stop - 1L)]
}

test_that("generic distribution CVLS routes use the checked block planner", {
  jksum_file <- locate_block_plan_source("jksum.c")
  planner_c_file <- locate_block_plan_source("jksum_block_plan.c")
  planner_h_file <- locate_block_plan_source("jksum_block_plan.h")
  skip_if(
    any(vapply(
      list(jksum_file, planner_c_file, planner_h_file),
      is.null,
      logical(1L)
    )),
    "complete package C sources unavailable in this test context"
  )

  jksum <- readLines(jksum_file, warn = FALSE)
  planner_c <- readLines(planner_c_file, warn = FALSE)
  planner_h <- readLines(planner_h_file, warn = FALSE)
  one_body <- extract_native_function(
    jksum,
    "^double np_kernel_estimate_distribution_ls_cv\\(",
    "^int np_kernel_estimate_con_distribution_categorical_leave_one_out_ls_cv\\("
  )
  two_body <- extract_native_function(
    jksum,
    "^int np_kernel_estimate_con_distribution_categorical_leave_one_out_ls_cv\\(",
    "^int np_kernel_estimate_con_density_categorical_leave_one_out_ls_cv\\("
  )

  expect_equal(sum(grepl(
    "np_jksum_distribution_one_block_plan_or_die(",
    one_body,
    fixed = TRUE
  )), 1L)
  expect_equal(sum(grepl(
    "np_jksum_distribution_two_block_plan_or_die(",
    two_body,
    fixed = TRUE
  )), 1L)
  expect_equal(sum(grepl(
    "np_jksum_memfac_cells(memfac, &Nm)",
    jksum,
    fixed = TRUE
  )), 0L)
  expect_equal(sum(grepl(
    "np_jksum_distribution_one_full_plan_safe(",
    one_body,
    fixed = TRUE
  )), 1L)
  expect_equal(sum(grepl(
    "np_jksum_block_plan_status_message(status)",
    jksum,
    fixed = TRUE
  )), 2L)

  obsolete <- c(
    "num_obs_eval_alloc*(num_obs_train_alloc+1)",
    "const int64_t wx0 = Nm/",
    "const int64_t wy0 =",
    "wy0 <= 0",
    "num_obs_train_alloc*num_obs_train_alloc*sizeof(double)",
    "num_obs_eval_alloc*num_obs_train_alloc*sizeof(double)"
  )
  for (pattern in obsolete) {
    expect_false(any(grepl(
      pattern,
      c(one_body, two_body),
      fixed = TRUE
    )))
  }

  expect_false(any(grepl(
    "malloc(num_obs_train_alloc*num_obs_wy_alloc*sizeof(double))",
    two_body,
    fixed = TRUE
  )))
  expect_false(any(grepl(
    "malloc(num_obs_wx_alloc*num_obs_train_alloc*sizeof(double))",
    two_body,
    fixed = TRUE
  )))
  expect_gte(sum(grepl(
    "np_jksum_malloc_array3_or_die(",
    two_body,
    fixed = TRUE
  )), 4L)

  module <- c(planner_c, planner_h)
  expect_true(any(grepl("np_size_add_checked(", module, fixed = TRUE)))
  expect_true(any(grepl("np_size_mul_checked(", module, fixed = TRUE)))
  expect_true(any(grepl(
    "np_size_array_bytes_checked(",
    module,
    fixed = TRUE
  )))
  expect_true(any(grepl(
    "NP_JKSUM_BLOCK_PLAN_CAPACITY",
    module,
    fixed = TRUE
  )))
  expect_false(any(grepl("malloc(", module, fixed = TRUE)))
  expect_false(any(grepl("Rprintf", module, fixed = TRUE)))
  expect_false(any(grepl("MPI_", module, fixed = TRUE)))
})
