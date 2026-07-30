locate_engine_sources <- function() {
  roots <- c(
    testthat::test_path("..", ".."),
    testthat::test_path("..", "..", ".."),
    Sys.getenv("R_PACKAGE_DIR", ""),
    Sys.getenv("R_PACKAGE_SOURCE", ""),
    getwd(),
    file.path(getwd(), "..")
  )
  roots <- unique(roots[nzchar(roots)])
  roots <- roots[file.exists(file.path(roots, "src", "kernelcv.c"))]
  if (!length(roots)) {
    return(NULL)
  }
  roots[[1L]]
}

test_that("legacy density translation unit and engines are absent", {
  root <- locate_engine_sources()
  skip_if(is.null(root), "package sources unavailable")

  expect_false(file.exists(file.path(root, "src", "kernele.c")))
  expect_true(file.exists(file.path(root, "src", "quantile.c")))

  native_lines <- c(
    readLines(file.path(root, "src", "kernelcv.c"), warn = FALSE),
    readLines(file.path(root, "src", "quantile.c"), warn = FALSE),
    readLines(file.path(root, "src", "headers.h"), warn = FALSE)
  )
  retired <- c(
    "kernel_estimate_density_categorical",
    "kernel_estimate_density_categorical_leave_one_out_cv",
    "kernel_estimate_distribution_categorical",
    "kernel_estimate_con_density_categorical",
    "kernel_estimate_con_density_categorical_leave_one_out_cv",
    "kernel_estimate_con_distribution_categorical",
    "kernel_estimate_con_distribution_categorical_leave_one_out",
    "kernel_estimate_con_distribution_categorical_leave_one_out_ccdf",
    "kernel_estimate_con_density_categorical_gradient",
    "kernel_estimate_con_density_categorical_gradient_categorical",
    "kernel_estimate_con_distribution_categorical_gradient",
    "kernel_estimate_con_distribution_categorical_gradient_categorical",
    "kernel_estimate_density_categorical_convolution_cv",
    "kernel_estimate_con_density_categorical_convolution_cv",
    "cv_func_density_categorical_ml",
    "cv_func_density_categorical_ls",
    "cv_func_con_density_categorical_ml",
    "cv_func_con_density_categorical_ls",
    "np_cv_func_con_density_categorical_ls",
    "cv_func_con_distribution_categorical_ccdf"
  )
  for (symbol in retired) {
    definition <- paste0(
      "^(int|double|void)[[:space:]]+",
      symbol,
      "[[:space:]]*\\("
    )
    expect_false(any(grepl(definition, native_lines)))
  }

  quantile <- readLines(file.path(root, "src", "quantile.c"), warn = FALSE)
  expect_equal(
    sum(grepl(
      "^int kernel_estimate_con_distribution_categorical_no_mpi\\(",
      quantile
    )),
    1L
  )
  expect_equal(
    sum(grepl("^int kernel_estimate_quantile\\(", quantile)),
    1L
  )
})

test_that("conditional-density CVLS wrapper has no dormant second engine", {
  root <- locate_engine_sources()
  skip_if(is.null(root), "package sources unavailable")

  lines <- readLines(file.path(root, "src", "jksum.c"), warn = FALSE)
  start <- grep(
    "^int np_kernel_estimate_con_density_categorical_leave_one_out_ls_cv\\(",
    lines
  )
  stop <- grep(
    "^static void np_lp_power2_moments_from_kernel_row\\(",
    lines
  )
  expect_length(start, 1L)
  expect_length(stop, 1L)
  expect_gt(stop, start)

  wrapper <- lines[start:(stop - 1L)]
  expect_lt(length(wrapper), 60L)
  expect_equal(
    sum(grepl(
      "return np_conditional_density_cvls_lp_stream(vector_scale_factor, cv);",
      wrapper,
      fixed = TRUE
    )),
    1L
  )
  expect_true(all(vapply(
    c("BW_FIXED", "BW_GEN_NN", "BW_ADAP_NN"),
    function(topology) any(grepl(topology, wrapper, fixed = TRUE)),
    logical(1)
  )))
  expect_false(any(grepl(
    "NP_GateOverrideCtx|np_gate_ctx_|alloc_matd|np_kernel_weighted_sum",
    wrapper
  )))
})
