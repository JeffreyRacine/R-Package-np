locate_beta_bandwidth_config_sources <- function() {
  roots <- unique(c(
    testthat::test_path("..", ".."),
    testthat::test_path("..", "..", ".."),
    Sys.getenv("R_PACKAGE_DIR", ""),
    Sys.getenv("R_PACKAGE_SOURCE", ""),
    getwd(),
    file.path(getwd(), "..")
  ))
  roots <- roots[nzchar(roots)]
  for (root in roots) {
    paths <- file.path(root, "src", c("beta_bandwidth.c", "kernelb.c"))
    if (all(file.exists(paths)))
      return(paths)
  }
  NULL
}

test_that("beta NN bandwidth preparation uses the state-free metric owner", {
  paths <- locate_beta_bandwidth_config_sources()
  skip_if(is.null(paths), "package sources unavailable")

  beta_source <- paste(readLines(paths[[1L]], warn = FALSE), collapse = "\n")
  kernel_source <- readLines(paths[[2L]], warn = FALSE)
  owner <- np_test_extract_c_function(
    kernel_source, "np_kernel_bandwidth_continuous_nn"
  )
  mean_owner <- np_test_extract_c_function(
    kernel_source, "kernel_bandwidth_mean_ctx"
  )

  owner_calls <- gregexpr(
    "np_kernel_bandwidth_continuous_nn(", beta_source, fixed = TRUE
  )[[1L]]
  expect_equal(sum(owner_calls > 0L), 2L)
  expect_false(grepl("kernel_bandwidth_mean(", beta_source, fixed = TRUE))

  for (token in c("int_LARGE_SF", "nconfac_extern", "ncatfac_extern")) {
    expect_false(grepl(token, beta_source, fixed = TRUE))
    expect_false(grepl(token, owner, fixed = TRUE))
  }
  expect_equal(
    sum(gregexpr(
      "np_kernel_bandwidth_continuous_nn_into_ctx(", mean_owner, fixed = TRUE
    )[[1L]] > 0L),
    2L
  )
})

test_that("interleaved beta and ordinary NN calls retain exact state", {
  old <- options(np.extendednn = FALSE)
  on.exit(options(old), add = TRUE)

  train <- data.frame(x = c(0.03, 0.12, 0.28, 0.51, 0.79, 0.96))
  evaluation <- data.frame(x = c(0.08, 0.34, 0.68, 0.91))
  beta_call <- function() {
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
  gaussian_call <- function() {
    npksum(
      bws = 2,
      txdat = train,
      exdat = evaluation,
      bwtype = "adaptive_nn",
      ckertype = "gaussian",
      return.kernel.weights = TRUE
    )
  }
  invalid_beta_call <- function() {
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

  beta_before <- beta_call()
  gaussian_before <- gaussian_call()
  expect_error(
    invalid_beta_call(),
    "invalid beta nearest-neighbor bandwidth"
  )
  gaussian_after <- gaussian_call()
  beta_after <- beta_call()

  expect_identical(gaussian_after$ksum, gaussian_before$ksum)
  expect_identical(gaussian_after$kw, gaussian_before$kw)
  expect_identical(beta_after$ksum, beta_before$ksum)
  expect_identical(beta_after$kw, beta_before$kw)
})
