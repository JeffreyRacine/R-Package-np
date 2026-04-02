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

test_that("linalg helpers clean up audited allocations before malloc errors", {
  src_file <- locate_native_source("linalg.c")
  skip_if(is.null(src_file), "source file src/linalg.c unavailable in this test context")

  src <- paste(readLines(src_file, warn = FALSE), collapse = "\n")

  expect_true(grepl("while \\(--i >= 0\\)\\s*free\\(\\*\\(\\(double \\*\\*\\)\\(&mat->matrix\\) \\+ i\\)\\);\\s*free\\(mat\\);\\s*error\\(\"mat: malloc error\\\\n\" \\);", src))
  expect_true(grepl("free\\(A\\);\\s*free\\(ipiv\\);\\s*error\\(\"mat_inv: malloc error\\\\n\"\\);", src))
  expect_true(grepl("free\\(Ac\\);\\s*free\\(Bc\\);\\s*free\\(ipiv\\);\\s*error\\(\"mat_solve: malloc error\\\\n\"\\);", src))
  expect_true(grepl("free\\(Ac\\);\\s*free\\(ipiv\\);\\s*error\\(\"mat_is_nonsingular: malloc error\\\\n\"\\);", src))
  expect_true(grepl("free\\(Ac\\);\\s*free\\(bc\\);\\s*free\\(ipiv\\);\\s*error\\(\"mat_inv00: malloc error\\\\n\"\\);", src))
})
