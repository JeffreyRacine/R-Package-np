locate_npcdist_resident_source <- function() {
  roots <- c(
    test_path("..", ".."),
    test_path("..", "..", ".."),
    Sys.getenv("R_PACKAGE_DIR", ""),
    Sys.getenv("R_PACKAGE_SOURCE", ""),
    getwd(),
    file.path(getwd(), "..")
  )
  roots <- unique(roots[nzchar(roots)])
  roots <- roots[file.exists(file.path(roots, "src", "jksum.c"))]
  if (!length(roots)) return(NULL)
  file.path(roots[[1L]], "src", "jksum.c")
}

test_that("all-large conditional-distribution CVLS keeps its objective resident", {
  path <- locate_npcdist_resident_source()
  skip_if(is.null(path), "package C source unavailable in this test context")
  source <- paste(readLines(path, warn = FALSE), collapse = "\n")
  start <- regexpr(
    "static int np_conditional_distribution_cvls_lp_all_large_stream",
    source,
    fixed = TRUE
  )[[1L]]
  stop <- regexpr(
    "static int np_conditional_density_cvml_lp_block_stream",
    source,
    fixed = TRUE
  )[[1L]]
  expect_gt(start, 0L)
  expect_gt(stop, start)
  body <- substr(source, start, stop - 1L)

  expect_match(body, "double cv_accumulator = 0.0;", fixed = TRUE)
  expect_match(body, "int cv_started = 0;", fixed = TRUE)
  expect_match(body, "cv_accumulator += tvd*tvd;", fixed = TRUE)
  expect_match(
    body,
    "if(np_distribution_cvls_finalize(\n       cv_accumulator,",
    fixed = TRUE
  )
  expect_match(
    body,
    "&cv_accumulator) != NP_DISTRIBUTION_CVLS_FINALIZE_OK)",
    fixed = TRUE
  )
  expect_match(
    body,
    "if(cv_started && status == 0)\n    *cv = cv_accumulator;",
    fixed = TRUE
  )
  expect_false(grepl("*cv += tvd*tvd;", body, fixed = TRUE))
})

test_that("empirical CDF-CV normalization has one checked finalizer", {
  path <- locate_npcdist_resident_source()
  skip_if(is.null(path), "package C source unavailable in this test context")
  source <- paste(readLines(path, warn = FALSE), collapse = "\n")

  expect_match(
    source,
    "static NPDistributionCvlsFinalizeStatus np_distribution_cvls_finalize(",
    fixed = TRUE
  )
  calls <- gregexpr(
    "np_distribution_cvls_finalize(", source, fixed = TRUE
  )[[1L]]
  expect_length(calls[calls > 0L], 23L)
  expect_false(grepl(
    "*cv /= np_distribution_cvls_pair_count(", source, fixed = TRUE
  ))
  expect_false(grepl(
    "*cv /= np_conditional_distribution_cvls_pair_count(", source,
    fixed = TRUE
  ))
  expect_false(grepl(
    "cv_accumulator /= np_conditional_distribution_cvls_pair_count(",
    source, fixed = TRUE
  ))
})
