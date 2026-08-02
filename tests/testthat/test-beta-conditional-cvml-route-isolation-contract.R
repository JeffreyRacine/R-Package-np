locate_beta_conditional_cvml_sources <- function() {
  roots <- unique(c(
    Sys.getenv("NP_SOURCE_ROOT", unset = ""),
    normalizePath(file.path(getwd(), "..", ".."), mustWork = FALSE),
    normalizePath(getwd(), mustWork = FALSE)
  ))
  roots <- roots[nzchar(roots)]
  for (root in roots) {
    if (file.exists(file.path(root, "src", "headers.h")) &&
        file.exists(file.path(root, "src", "jksum.c")) &&
        file.exists(file.path(root, "src", "kernelcv.c")))
      return(root)
  }
  NULL
}

beta_conditional_count_fixed <- function(text, pattern) {
  hits <- gregexpr(pattern, text, fixed = TRUE)[[1L]]
  sum(hits > 0L)
}

test_that("conditional CVML route metadata has one typed X/Y contract", {
  root <- locate_beta_conditional_cvml_sources()
  skip_if(is.null(root), "package sources unavailable")
  headers <- paste(
    readLines(file.path(root, "src", "headers.h"), warn = FALSE),
    collapse = "\n"
  )

  expect_match(
    headers,
    paste0(
      "typedef struct {\n",
      "  const NPContinuousKernelRoute *x_route;\n",
      "  NPContinuousKernelDerivativeDiagnostics *x_diagnostics;\n",
      "  const NPContinuousKernelRoute *y_route;\n",
      "  NPContinuousKernelDerivativeDiagnostics *y_diagnostics;\n",
      "  int categorical_compress;\n",
      "} NPConditionalKernelExecutionContext;"
    ),
    fixed = TRUE
  )
  expect_match(
    headers,
    "np_kernel_estimate_con_density_categorical_leave_one_out_cv_ctx(",
    fixed = TRUE
  )
  expect_match(
    headers,
    "const NPConditionalKernelExecutionContext *execution_context,",
    fixed = TRUE
  )
})

test_that("conditional CVML route sibling is dormant beside a literal incumbent", {
  root <- locate_beta_conditional_cvml_sources()
  skip_if(is.null(root), "package sources unavailable")
  engine <- paste(
    readLines(file.path(root, "src", "jksum.c"), warn = FALSE),
    collapse = "\n"
  )
  ingress <- paste(
    readLines(file.path(root, "src", "kernelcv.c"), warn = FALSE),
    collapse = "\n"
  )

  incumbent_start <- regexpr(
    "int np_kernel_estimate_con_density_categorical_leave_one_out_cv(",
    engine,
    fixed = TRUE
  )[[1L]]
  incumbent_end <- regexpr(
    "static int np_conditional_categorical_profile_fit(",
    engine,
    fixed = TRUE
  )[[1L]]
  sibling_start <- regexpr(
    "int np_kernel_estimate_con_density_categorical_leave_one_out_cv_ctx(",
    engine,
    fixed = TRUE
  )[[1L]]
  expect_gt(incumbent_start, 0L)
  expect_gt(incumbent_end, incumbent_start)
  expect_gt(sibling_start, incumbent_end)

  incumbent <- substr(engine, incumbent_start, incumbent_end - 1L)
  expect_false(grepl("NPConditionalKernelExecutionContext", incumbent,
                     fixed = TRUE))
  expect_false(grepl("_cv_ctx(", ingress, fixed = TRUE))
  expect_identical(
    beta_conditional_count_fixed(
      ingress,
      "np_kernel_estimate_con_density_categorical_leave_one_out_cv("
    ),
    1L
  )
  expect_match(
    ingress,
    "np_beta_conditional_density_bw_objective(vector_scale_factor, 0)",
    fixed = TRUE
  )
  expect_false(grepl(
    "np_conditional_density_bw_categorical_compress_extern",
    ingress,
    fixed = TRUE
  ))
})

test_that("conditional CVML route sibling delegates null and fails closed", {
  root <- locate_beta_conditional_cvml_sources()
  skip_if(is.null(root), "package sources unavailable")
  engine <- paste(
    readLines(file.path(root, "src", "jksum.c"), warn = FALSE),
    collapse = "\n"
  )
  sibling_start <- regexpr(
    "int np_kernel_estimate_con_density_categorical_leave_one_out_cv_ctx(",
    engine,
    fixed = TRUE
  )[[1L]]
  expect_gt(sibling_start, 0L)
  sibling <- substr(engine, sibling_start, nchar(engine))

  expect_match(
    sibling,
    paste0(
      "if(execution_context == NULL)\n",
      "    return np_kernel_estimate_con_density_categorical_leave_one_out_cv("
    ),
    fixed = TRUE
  )
  expect_match(
    sibling,
    "np_conditional_kernel_execution_context_valid(",
    fixed = TRUE
  )
  expect_match(
    sibling,
    "conditional density CVML kernel route has an invalid layout",
    fixed = TRUE
  )
  expect_match(sibling, "return 1;", fixed = TRUE)
  expect_false(grepl("kernel_weighted_sum_np(", sibling, fixed = TRUE))
  expect_false(grepl("malloc(", sibling, fixed = TRUE))
  expect_false(grepl("calloc(", sibling, fixed = TRUE))
  expect_false(grepl("realloc(", sibling, fixed = TRUE))

  expect_match(
    engine,
    "np_continuous_kernel_route_validate(",
    fixed = TRUE
  )
  expect_match(
    engine,
    "np_continuous_kernel_route_has_beta(",
    fixed = TRUE
  )
  expect_match(
    engine,
    "context->categorical_compress != 0",
    fixed = TRUE
  )
  expect_match(
    engine,
    "context->x_route->segment_count",
    fixed = TRUE
  )
})
