library(np)

locate_native_source <- function(filename) {
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
  if (length(hits) == 0L) {
    return(NULL)
  }
  hits[[1L]]
}

test_that("tree allocation hotspots use overflow-checked helpers", {
  src_file <- locate_native_source("tree.c")
  skip_if(is.null(src_file), "source file src/tree.c unavailable in this test context")

  src <- paste(readLines(src_file, warn = FALSE), collapse = "\n")

  expect_true(grepl("np_tree_size_mul_or_die", src, fixed = TRUE))
  expect_true(grepl("realloc(ptr, np_tree_size_mul_or_die(", src, fixed = TRUE))
  expect_true(grepl("np_tree_malloc_array_or_die((size_t)numnode, sizeof(KDN), \"build_kdtree kdn\")", src, fixed = TRUE))
  expect_true(grepl("np_tree_size_mul3_or_die", src, fixed = TRUE))
  expect_false(grepl("kdx->kdn = \\(KDN \\*\\)malloc\\(numnode\\*sizeof\\(KDN\\)\\);", src))
  expect_false(grepl("kdx->bb = \\(double \\*\\)malloc\\(numnode\\*2\\*ndim\\*sizeof\\(double\\)\\);", src))
})

test_that("audited jksum allocation hotspots use checked helpers", {
  src_file <- locate_native_source("jksum.c")
  skip_if(is.null(src_file), "source file src/jksum.c unavailable in this test context")

  src <- paste(readLines(src_file, warn = FALSE), collapse = "\n")

  expect_true(grepl("np_jksum_malloc_array_or_die", src, fixed = TRUE))
  expect_true(grepl("np_jksum_malloc_array3_or_die", src, fixed = TRUE))
  expect_true(grepl("np_kernel_estimate_density_categorical_leave_one_out_cv kwx", src, fixed = TRUE))
  expect_true(grepl("np_kernel_estimate_density_categorical_convolution_cv operator", src, fixed = TRUE))
  expect_true(grepl("np_kernel_estimate_con_dens_dist_categorical permn", src, fixed = TRUE))
  expect_false(grepl("permn = \\(double \\*\\)malloc\\(num_X\\*num_obs_eval_alloc\\*sizeof\\(double\\)\\);", src))
  expect_false(grepl("permd = \\(double \\*\\)malloc\\(num_X\\*num_obs_eval_alloc\\*sizeof\\(double\\)\\);", src))
})
