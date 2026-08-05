locate_beta_facade_sources <- function() {
  roots <- unique(c(
    testthat::test_path("..", ".."),
    testthat::test_path("..", "..", ".."),
    Sys.getenv("R_PACKAGE_DIR", ""),
    Sys.getenv("R_PACKAGE_SOURCE", ""),
    getwd(),
    file.path(getwd(), "..")
  ))
  roots <- roots[nzchar(roots)]
  roots <- roots[file.exists(file.path(roots, "src", "beta_bandwidth.h"))]
  if (!length(roots))
    return(NULL)
  roots[[1L]]
}

test_that("the dead beta kernel-sum facade is absent", {
  root <- locate_beta_facade_sources()
  skip_if(is.null(root), "package sources unavailable")

  source_dir <- file.path(root, "src")
  expect_false(file.exists(file.path(source_dir, "beta_kernelsum.c")))
  expect_false(file.exists(file.path(source_dir, "beta_kernelsum.h")))

  bandwidth_contract <- paste(
    readLines(file.path(source_dir, "beta_bandwidth.h"), warn = FALSE),
    collapse = "\n"
  )
  expect_match(bandwidth_contract, "NP_BETA_BANDWIDTH_FIXED = 0", fixed = TRUE)
  expect_match(
    bandwidth_contract,
    "NP_BETA_BANDWIDTH_GENERALIZED_NN = 1",
    fixed = TRUE
  )
  expect_match(
    bandwidth_contract,
    "NP_BETA_BANDWIDTH_ADAPTIVE_NN = 2",
    fixed = TRUE
  )

  native_files <- list.files(
    source_dir,
    pattern = "[.][ch]$",
    full.names = TRUE
  )
  native_source <- paste(
    unlist(lapply(native_files, readLines, warn = FALSE), use.names = FALSE),
    collapse = "\n"
  )
  expect_false(grepl("#include \"beta_kernelsum.h\"", native_source,
                     fixed = TRUE))
  expect_false(grepl("np_beta_kernelsum(", native_source, fixed = TRUE))

  retired_kernel_roots <- c(
    "np_beta_pdf_scale(",
    "np_beta_pdf_order2(",
    "np_beta_pdf_order(",
    "np_beta_log_pdf_order2(",
    "np_beta_derivative_regular_value(",
    "np_beta_derivative_public_value(",
    "np_beta_overlap_order2(",
    "np_beta_overlap_order("
  )
  for (root_name in retired_kernel_roots) {
    expect_false(
      grepl(root_name, native_source, fixed = TRUE),
      info = root_name
    )
  }
  expect_match(
    native_source,
    "np_beta_log_abs_pdf_order(",
    fixed = TRUE
  )
  expect_match(
    native_source,
    "np_beta_log_abs_overlap_order(",
    fixed = TRUE
  )
  expect_equal(
    lengths(regmatches(
      native_source,
      gregexpr("NP_BETA_BANDWIDTH_FIXED = 0", native_source, fixed = TRUE)
    )),
    1L
  )
})
