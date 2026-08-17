locate_kernel_bandwidth_source <- function() {
  roots <- unique(c(
    testthat::test_path("..", ".."),
    testthat::test_path("..", "..", ".."),
    Sys.getenv("R_PACKAGE_DIR", ""),
    Sys.getenv("R_PACKAGE_SOURCE", ""),
    getwd(),
    file.path(getwd(), "..")
  ))
  roots <- roots[nzchar(roots)]
  paths <- file.path(roots, "src", "kernelb.c")
  paths <- paths[file.exists(paths)]
  if (!length(paths))
    return(NULL)
  paths[[1L]]
}

test_that("shared bandwidth owners clean ordinary NN failures", {
  path <- locate_kernel_bandwidth_source()
  skip_if(is.null(path), "package sources unavailable")
  source <- readLines(path, warn = FALSE)

  regions <- list(
    kernel_bandwidth = list(
      code = npRmpi_test_extract_c_function(source, "kernel_bandwidth_ctx"),
      cleanup_gotos = 8L
    ),
    kernel_bandwidth_mean = list(
      code = npRmpi_test_extract_c_function(
        source, "kernel_bandwidth_mean_ctx"
      ),
      cleanup_gotos = 2L
    )
  )

  for (contract in regions) {
    region <- contract$code
    allocation <- regexpr("nn_distance = alloc_vecd", region, fixed = TRUE)
    cleanup <- regexpr("\ncleanup:", region, fixed = TRUE)
    expect_gt(allocation, 0L)
    expect_gt(cleanup, allocation)
    owned_region <- substr(region, allocation, cleanup)
    expect_false(grepl("return(1);", owned_region, fixed = TRUE))
    gotos <- gregexpr("goto cleanup;", region, fixed = TRUE)[[1L]]
    expect_equal(sum(gotos > 0L), contract$cleanup_gotos)
    expect_match(region, "free(nn_distance);", fixed = TRUE)
    expect_match(region, "free(vec_sdev_x);", fixed = TRUE)
    expect_match(region, "free(vec_sdev_y);", fixed = TRUE)
    expect_match(region, "return(status);", fixed = TRUE)
  }

  owner <- npRmpi_test_extract_c_function(
    source, "np_kernel_bandwidth_continuous_nn"
  )
  expect_match(owner, "nn_distance = alloc_vecd", fixed = TRUE)
  expect_match(owner, "free(nn_distance);", fixed = TRUE)
  expect_match(owner, "return status;", fixed = TRUE)
})

test_that("failed beta NN preparation is followed by an exact valid call", {
  active <- isTRUE(try(
    getFromNamespace(".npRmpi_has_active_slave_pool", "npRmpi")(),
    silent = TRUE
  ))
  skip_if_not(active, "runtime contract requires a test-owned MPI pool")

  old <- options(np.extendednn = FALSE)
  on.exit(options(old), add = TRUE)

  train <- data.frame(x = c(0.03, 0.12, 0.28, 0.51, 0.79, 0.96))
  evaluation <- data.frame(x = c(0.08, 0.34, 0.68, 0.91))
  valid_call <- function() {
    npksum(
      bws = 2,
      txdat = train,
      exdat = evaluation,
      bwtype = "adaptive_nn",
      ckertype = "beta",
      ckerorder = 8,
      ckerbound = "fixed",
      ckerlb = 0,
      ckerub = 1,
      return.kernel.weights = TRUE
    )
  }
  invalid_call <- function() {
    npksum(
      bws = 7,
      txdat = train,
      exdat = evaluation,
      bwtype = "adaptive_nn",
      ckertype = "beta",
      ckerbound = "fixed",
      ckerlb = 0,
      ckerub = 1
    )
  }

  before <- valid_call()
  for (i in seq_len(25L))
    expect_error(invalid_call(), "invalid beta nearest-neighbor bandwidth")
  after <- valid_call()

  expect_identical(after$ksum, before$ksum)
  expect_identical(after$kw, before$kw)
})
