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
  if (!length(hits)) {
    return(NULL)
  }
  hits[[1L]]
}

test_that("canonical regression objective sizes are checked before allocation", {
  path <- locate_native_source("jksum.c")
  skip_if(is.null(path), "source file src/jksum.c unavailable")
  source <- paste(readLines(path, warn = FALSE), collapse = "\n")

  expect_true(grepl('#include "np_native_safety.h"', source, fixed = TRUE))
  expect_true(grepl("np_size_add_checked((size_t)nrc1", source, fixed = TRUE))
  expect_true(grepl("np_size_mul_checked(nrc2_size, nrc2_size", source, fixed = TRUE))
  expect_true(grepl("np_native_malloc_array((void **)&mean", source, fixed = TRUE))
  expect_true(grepl("np_native_malloc_array((void **)&lc_Y[1]", source, fixed = TRUE))
  expect_false(grepl("malloc(2*num_obs_eval_alloc*sizeof(double))", source, fixed = TRUE))
  expect_false(grepl("const int nrcc22 = nrc2*nrc2", source, fixed = TRUE))
})

test_that("Gaussian MPI partition arithmetic is checked before narrowing", {
  path <- locate_native_source("jksum_gaussian_density.c")
  skip_if(is.null(path), "source file src/jksum_gaussian_density.c unavailable")
  source <- paste(readLines(path, warn = FALSE), collapse = "\n")

  expect_true(grepl("np_int_ceil_div_nonnegative", source, fixed = TRUE))
  expect_true(grepl("np_size_mul_checked((size_t)stride", source, fixed = TRUE))
  expect_true(grepl("np_size_to_int_checked(eval_start_size", source, fixed = TRUE))
  expect_false(
    grepl("(num_obs + iNum_Processors - 1)/iNum_Processors", source, fixed = TRUE)
  )
})
