library(npRmpi)

locate_r_source <- function(filename) {
  candidates <- c(
    test_path("..", "..", "R", filename),
    test_path("..", "..", "..", "R", filename),
    file.path(Sys.getenv("R_PACKAGE_DIR", ""), "R", filename),
    file.path(Sys.getenv("R_PACKAGE_SOURCE", ""), "R", filename),
    file.path(getwd(), "R", filename),
    file.path(getwd(), "..", "R", filename)
  )
  candidates <- unique(candidates[nzchar(candidates)])
  hits <- candidates[file.exists(candidates)]
  if (length(hits) == 0L)
    return("")
  hits[[1L]]
}

test_that("timing producers avoid user CPU summaries for touched families", {
  files <- c(
    "np.density.bw.R",
    "np.distribution.bw.R",
    "np.condensity.bw.R",
    "np.condistribution.bw.R",
    "np.plregression.bw.R",
    "np.singleindex.bw.R",
    "np.smoothcoef.bw.R"
  )

  paths <- vapply(files, locate_r_source, character(1), USE.NAMES = TRUE)
  if (any(!nzchar(paths)))
    skip("R source files unavailable in this test context")

  for (nm in names(paths)) {
    lines <- readLines(paths[[nm]], warn = FALSE)
    expect_false(any(grepl("\\}\\)\\[1\\]", lines)), info = nm)
  }
})

test_that("reference regression bandwidth timing remains full-stage elapsed", {
  src_file <- locate_r_source("np.regression.bw.R")
  skip_if(!nzchar(src_file) || !file.exists(src_file),
          "R source file np.regression.bw.R unavailable in this test context")

  lines <- readLines(src_file, warn = FALSE)
  if (!any(grepl("npregbw", lines, fixed = TRUE)))
    skip("R source file np.regression.bw.R unavailable in this test context")
  expect_true(any(grepl("elapsed.start <- proc.time\\(\\)\\[3\\]", lines)), info = "np.regression.bw.R")
  expect_true(any(grepl("tbw\\$total.time <- proc.time\\(\\)\\[3\\] - elapsed.start", lines)),
              info = "np.regression.bw.R")
})
