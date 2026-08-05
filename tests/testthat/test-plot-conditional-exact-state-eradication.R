locate_conditional_bootstrap_sources <- function() {
  roots <- unique(c(
    testthat::test_path("..", ".."),
    testthat::test_path("..", "..", ".."),
    Sys.getenv("R_PACKAGE_DIR", ""),
    Sys.getenv("R_PACKAGE_SOURCE", ""),
    getwd(),
    file.path(getwd(), "..")
  ))
  roots <- roots[nzchar(roots)]
  roots <- roots[file.exists(file.path(roots, "R", "np.plot.helpers.R"))]
  if (!length(roots))
    return(NULL)
  roots[[1L]]
}

test_that("conditional bootstrap exact-state sidecars are absent", {
  root <- locate_conditional_bootstrap_sources()
  skip_if(is.null(root), "package sources unavailable")

  helpers <- paste(
    readLines(file.path(root, "R", "np.plot.helpers.R"), warn = FALSE),
    collapse = "\n"
  )
  expect_false(grepl(".np_ksum_exact_state_build", helpers, fixed = TRUE))
  expect_false(grepl(".np_ksum_eval_exact_state", helpers, fixed = TRUE))
  expect_false(grepl(
    ".np_ksum_conditional_eval_exact_boot_active",
    helpers,
    fixed = TRUE
  ))
  expect_match(
    helpers,
    ".np_inid_boot_from_ksum_conditional_exact <- function",
    fixed = TRUE
  )
  expect_match(
    helpers,
    ".np_inid_boot_from_conditional_localpoly_fixed",
    fixed = TRUE
  )
})
