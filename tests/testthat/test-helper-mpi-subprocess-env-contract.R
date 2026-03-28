test_that("npRmpi_subprocess_env does not treat the source tree parent as a library root", {
  env <- npRmpi_subprocess_env()
  skip_if(is.null(env), "subprocess library setup unavailable")

  rlibs <- env[grepl("^R_LIBS=", env)]
  expect_identical(length(rlibs), 1L)

  lib.paths <- strsplit(sub("^R_LIBS=", "", rlibs), .Platform$path.sep, fixed = TRUE)[[1L]]
  pkg.root <- normalizePath(test_path("..", ".."), mustWork = TRUE)
  bad.root <- normalizePath(dirname(pkg.root), mustWork = TRUE)

  expect_false(any(vapply(
    lib.paths,
    function(path) identical(normalizePath(path, mustWork = FALSE), bad.root),
    logical(1)
  )))
})
