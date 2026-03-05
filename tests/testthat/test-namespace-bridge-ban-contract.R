test_that("npRmpi repo contains no serial namespace-bridge tokens", {
  root <- normalizePath(testthat::test_path("..", ".."), mustWork = TRUE)
  dirs <- file.path(root, c("R", "tests", "demo", "vignettes", "man", "issue_notes"))
  dirs <- dirs[dir.exists(dirs)]
  skip_if(length(dirs) == 0L, "source tree unavailable in installed test context")

  files <- unlist(lapply(dirs, function(d) list.files(d, recursive = TRUE, full.names = TRUE)))
  files <- files[file.exists(files) & !dir.exists(files)]

  # Text-like files only.
  text.ext <- c(
    ".R", ".Rmd", ".Rnw", ".Rd", ".md", ".txt", ".sh",
    ".c", ".h", ".f", ".f90", ".for"
  )
  keep <- vapply(files, function(f) {
    ext <- tools::file_ext(f)
    if (!nzchar(ext)) return(FALSE)
    paste0(".", tolower(ext)) %in% text.ext
  }, logical(1))
  files <- files[keep]

  bridge2 <- paste0("np", ":", ":")
  bridge3 <- paste0("np", ":", ":", ":")
  offenders <- character()

  for (f in files) {
    raw <- readLines(f, warn = FALSE)
    hit <- which(grepl(bridge2, raw, fixed = TRUE) | grepl(bridge3, raw, fixed = TRUE))
    if (length(hit)) {
      offenders <- c(
        offenders,
        sprintf("%s:%d: %s", sub(paste0("^", root, "/"), "", f), hit, trimws(raw[hit]))
      )
    }
  }

  expect_equal(length(offenders), 0, info = paste(offenders, collapse = "\n"))
})
