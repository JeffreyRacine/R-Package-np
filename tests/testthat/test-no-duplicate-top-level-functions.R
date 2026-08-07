top_level_function_bindings <- function(r_dir) {
  files <- sort(list.files(r_dir, pattern = "\\.[Rr]$", full.names = TRUE))
  bindings <- lapply(files, function(file) {
    exprs <- parse(file, keep.source = FALSE)
    Filter(nzchar, vapply(exprs, function(expr) {
      assignment <- is.call(expr) &&
        as.character(expr[[1L]]) %in% c("<-", "=") &&
        length(expr) == 3L && is.symbol(expr[[2L]]) &&
        is.call(expr[[3L]]) && identical(expr[[3L]][[1L]], as.name("function"))
      if (assignment) as.character(expr[[2L]]) else ""
    }, character(1L)))
  })
  unlist(bindings, use.names = FALSE)
}

test_that("package R sources have one top-level owner per function name", {
  r_dir <- normalizePath(testthat::test_path("..", "..", "R"),
                         mustWork = FALSE)
  skip_if_not(dir.exists(r_dir), "package source tree is unavailable")
  bindings <- top_level_function_bindings(r_dir)
  expect_gt(length(bindings), 0L)
  expect_identical(anyDuplicated(bindings), 0L,
                   info = paste(unique(bindings[duplicated(bindings)]), collapse = ", "))
})
