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

test_that("mat_vec allocation helpers free partial state before error", {
  src_file <- locate_native_source("mat_vec.c")
  skip_if(is.null(src_file), "source file src/mat_vec.c unavailable in this test context")

  src <- paste(readLines(src_file, warn = FALSE), collapse = "\n")

  expect_true(grepl("free(m);", src, fixed = TRUE))
  expect_true(grepl("free_mat(m, i);", src, fixed = TRUE))
})

test_that("the retired legacy linear-algebra source does not return", {
  marker <- locate_native_source("mat_vec.c")
  skip_if(is.null(marker), "package native source unavailable in this test context")
  src_dir <- dirname(marker)

  expect_false(file.exists(file.path(src_dir, paste0("lin", "alg.c"))))
  expect_false(file.exists(file.path(src_dir, paste0("lin", "alg.h"))))
})
