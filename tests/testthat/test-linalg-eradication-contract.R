library(npRmpi)

locate_native_source_for_linalg_contract <- function(filename) {
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
  if (!length(hits))
    return(NULL)
  hits[[1L]]
}

test_that("legacy linear-algebra files, types, and calls are absent", {
  marker <- locate_native_source_for_linalg_contract("jksum.c")
  skip_if(is.null(marker), "package native source unavailable in this test context")
  src_dir <- dirname(marker)

  expect_false(file.exists(file.path(src_dir, "linalg.c")))
  expect_false(file.exists(file.path(src_dir, "linalg.h")))

  native_files <- list.files(
    src_dir,
    pattern = "[.](c|h)$",
    full.names = TRUE
  )
  source <- paste(
    unlist(lapply(native_files, readLines, warn = FALSE), use.names = FALSE),
    collapse = "\n"
  )
  forbidden <- c(
    '#include "linalg.h"',
    "typedef double **MATRIX",
    "mat_creat(",
    "mat_fill(",
    "mat_free(",
    "mat_inv(",
    "mat_solve(",
    "mat_is_nonsingular(",
    "mat_inv00(",
    "isFiniteMatrix("
  )
  for (token in forbidden)
    expect_false(grepl(token, source, fixed = TRUE), info = token)
})

test_that("the final scalar and contiguous replacement owners remain explicit", {
  src_file <- locate_native_source_for_linalg_contract("jksum.c")
  skip_if(is.null(src_file), "source file src/jksum.c unavailable in this test context")
  src <- paste(readLines(src_file, warn = FALSE), collapse = "\n")

  expect_true(grepl("xtxinv = 1.0/xtx;", src, fixed = TRUE))
  expect_true(grepl("double *moment = NULL, *temp = NULL;", src, fixed = TRUE))
  expect_true(grepl("double *quad_mat = NULL;", src, fixed = TRUE))
  expect_true(grepl(
    "moment[a*ctx->nterms+b] += za*cross_terms[b];",
    src,
    fixed = TRUE
  ))
})
