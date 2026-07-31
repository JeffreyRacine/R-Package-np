locate_objective_cache_source <- function(filename) {
  candidates <- c(
    test_path("..", "..", "src", filename),
    test_path("..", "..", "..", "src", filename),
    file.path(Sys.getenv("R_PACKAGE_DIR", ""), "src", filename),
    file.path(Sys.getenv("R_PACKAGE_SOURCE", ""), "src", filename),
    file.path(getwd(), "src", filename),
    file.path(getwd(), "..", "src", filename)
  )
  candidates <- unique(candidates[nzchar(candidates)])
  hits <- candidates[file.exists(candidates)]
  if (!length(hits)) NULL else hits[[1L]]
}

test_that("native objective caches have one checked peak-memory contract", {
  source_file <- locate_objective_cache_source("np.c")
  header_file <- locate_objective_cache_source("np_native_safety.h")
  skip_if(is.null(source_file), "source file src/np.c unavailable")
  skip_if(is.null(header_file), "source file src/np_native_safety.h unavailable")
  source <- paste(readLines(source_file, warn = FALSE), collapse = "\n")
  header <- paste(readLines(header_file, warn = FALSE), collapse = "\n")

  expect_match(source, '#include "np_native_safety.h"', fixed = TRUE)
  expect_match(source, "#define NP_BWM_CACHE_MAX_BYTES", fixed = TRUE)
  expect_match(source, "(size_t)64U*(size_t)1024U*(size_t)1024U", fixed = TRUE)
  expect_match(header, "np_native_cache_table_bytes(", fixed = TRUE)
  expect_match(header, "np_native_cache_rehash_peak_bytes(", fixed = TRUE)
  expect_match(header, "np_native_cache_growth_size(", fixed = TRUE)
  expect_match(source, "bwm_nn_cache_growth_size", fixed = TRUE)
  expect_match(source, "bwm_objective_cache_growth_size", fixed = TRUE)
  expect_match(source, "peak_bytes <= NP_BWM_CACHE_MAX_BYTES", fixed = TRUE)
  expect_match(source, "bwm_cache_key_length_checked(", fixed = TRUE)
  expect_false(grepl(
    "capacity * (size_t)bwm_",
    source,
    fixed = TRUE
  ))
  expect_false(grepl("_cache_capacity * 2", source, fixed = TRUE))
})

test_that("native cache growth is transactional and rank-safe", {
  source_file <- locate_objective_cache_source("np.c")
  skip_if(is.null(source_file), "source file src/np.c unavailable")
  source <- paste(readLines(source_file, warn = FALSE), collapse = "\n")

  expect_false(grepl("bwm_cache_collective_available(", source, fixed = TRUE))
  expect_match(
    source,
    "bwm_nn_cache_requested && comm[1] != MPI_COMM_NULL",
    fixed = TRUE
  )
  expect_match(source, "&global_hit,", fixed = TRUE)
  expect_match(
    source,
    "1,\n                  MPI_INT,\n                  MPI_MIN",
    fixed = TRUE
  )
  expect_match(source, "bwm_cache_sync_failure_stats(", fixed = TRUE)
  expect_match(
    source,
    "failed,\n                  2,\n                  MPI_INT,\n                  MPI_MAX",
    fixed = TRUE
  )
  expect_match(source, "bwm_nn_cache_insert_enabled = 0", fixed = TRUE)
  expect_match(source, "bwm_objective_cache_insert_enabled = 0", fixed = TRUE)
  expect_match(source, "if (!bwm_nn_cache_insert_enabled)", fixed = TRUE)
  expect_match(source, "if (!bwm_objective_cache_insert_enabled)", fixed = TRUE)
  expect_match(source, "int *new_keys = NULL", fixed = TRUE)
  expect_match(source, "double *new_keys = NULL", fixed = TRUE)
  expect_false(grepl(
    "bwm_nn_cache_keys = NULL;\n  bwm_nn_cache_values = NULL;\n  bwm_nn_cache_used = NULL;\n  bwm_nn_cache_capacity = 0;\n  bwm_nn_cache_size = 0;\n  if (!bwm_nn_cache_alloc_table",
    source,
    fixed = TRUE
  ))
  expect_false(grepl(
    "bwm_objective_cache_keys = NULL;\n  bwm_objective_cache_values = NULL;\n  bwm_objective_cache_used = NULL;\n  bwm_objective_cache_capacity = 0;\n  bwm_objective_cache_size = 0;\n  if (!bwm_objective_cache_alloc_table",
    source,
    fixed = TRUE
  ))
})

test_that("wide nearest-neighbor cache keys use lifecycle-owned scratch", {
  source_file <- locate_objective_cache_source("np.c")
  skip_if(is.null(source_file), "source file src/np.c unavailable")
  source <- paste(readLines(source_file, warn = FALSE), collapse = "\n")

  expect_match(source, "bwm_nn_cache_dynamic_key", fixed = TRUE)
  expect_match(source, "np_native_malloc_array((void **)&dynamic_key", fixed = TRUE)
  expect_match(
    source,
    "(bwm_nn_cache_key_len > 16) ? bwm_nn_cache_dynamic_key : stack_key",
    fixed = TRUE
  )
  expect_false(grepl(
    "malloc((size_t)bwm_nn_cache_key_len * sizeof(int))",
    source,
    fixed = TRUE
  ))
})

test_that("np.objective.cache remains a strict logical option", {
  old <- options(np.objective.cache = "yes")
  on.exit(options(old), add = TRUE)
  expect_error(
    npRmpi:::npObjectiveCacheEnabled(),
    "must be TRUE or FALSE",
    fixed = TRUE
  )
})
