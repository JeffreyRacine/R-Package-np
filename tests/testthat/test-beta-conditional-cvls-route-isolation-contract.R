locate_beta_conditional_cvls_sources <- function() {
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

beta_conditional_cvls_count_fixed <- function(text, pattern) {
  hits <- gregexpr(pattern, text, fixed = TRUE)[[1L]]
  sum(hits > 0L)
}

test_that("conditional CVLS route sibling is the sole beta ingress", {
  root <- locate_beta_conditional_cvls_sources()
  skip_if(is.null(root), "package sources unavailable")
  headers <- paste(
    readLines(file.path(root, "src", "headers.h"), warn = FALSE),
    collapse = "\n"
  )
  engine <- paste(
    readLines(file.path(root, "src", "jksum.c"), warn = FALSE),
    collapse = "\n"
  )
  ingress <- paste(
    readLines(file.path(root, "src", "kernelcv.c"), warn = FALSE),
    collapse = "\n"
  )
  native.files <- list.files(
    file.path(root, "src"), pattern = "\\.[ch]$", full.names = TRUE
  )
  native <- paste(vapply(native.files, function(path) {
    paste(readLines(path, warn = FALSE), collapse = "\n")
  }, character(1)), collapse = "\n")

  expect_match(
    headers,
    "int np_conditional_density_cvls_lp_stream_ctx(",
    fixed = TRUE
  )
  expect_match(
    headers,
    "int np_conditional_distribution_cvls_lp_stream_ctx(",
    fixed = TRUE
  )
  expect_match(
    headers,
    "const NPConditionalKernelExecutionContext *execution_context,",
    fixed = TRUE
  )
  expect_identical(
    beta_conditional_cvls_count_fixed(
      engine,
      "int np_conditional_density_cvls_lp_stream_ctx("
    ),
    1L
  )
  expect_identical(
    beta_conditional_cvls_count_fixed(
      engine,
      "int np_conditional_distribution_cvls_lp_stream_ctx("
    ),
    1L
  )
  expect_identical(
    beta_conditional_cvls_count_fixed(
      ingress,
      "np_conditional_distribution_cvls_lp_stream_ctx("
    ),
    1L
  )
  expect_identical(
    beta_conditional_cvls_count_fixed(
      ingress,
      "np_conditional_density_cvls_lp_stream_ctx("
    ),
    1L
  )
  expect_match(
    ingress,
    "np_beta_conditional_density_bw_objective_ls_ctx(",
    fixed = TRUE
  )
  expect_match(
    ingress,
    "np_beta_conditional_distribution_bw_objective_ls_ctx(",
    fixed = TRUE
  )
  expect_false(grepl("np_beta_objective_conditional_density_ls(", ingress,
                     fixed = TRUE))
  expect_false(grepl(
    "np_beta_objective_conditional_density_ls(",
    native,
    fixed = TRUE
  ))
  expect_false(grepl(
    "np_beta_objective_conditional_density_ls_unbounded_legacy_y(",
    native,
    fixed = TRUE
  ))
  expect_false(grepl(
    "np_beta_objective_conditional_distribution_ls(",
    ingress,
    fixed = TRUE
  ))
  expect_false(file.exists(file.path(root, "src", "beta_objectives.c")))
  expect_false(file.exists(file.path(root, "src", "beta_objectives.h")))
  expect_false(grepl("np_beta_objective_", native, fixed = TRUE))
})

test_that("conditional CVLS route sibling delegates null and owns routed adapter", {
  root <- locate_beta_conditional_cvls_sources()
  skip_if(is.null(root), "package sources unavailable")
  engine <- paste(
    readLines(file.path(root, "src", "jksum.c"), warn = FALSE),
    collapse = "\n"
  )
  sibling_start <- regexpr(
    "int np_conditional_density_cvls_lp_stream_ctx(",
    engine,
    fixed = TRUE
  )[[1L]]
  expect_gt(sibling_start, 0L)
  sibling <- substr(engine, sibling_start, nchar(engine))

  expect_match(
    sibling,
    paste0(
      "if(execution_context == NULL)\n",
      "    return np_conditional_density_cvls_lp_stream("
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
    "conditional density CVLS kernel route has an invalid layout",
    fixed = TRUE
  )
  expect_match(sibling, "NPConditionalCVLSRowProvider provider;", fixed = TRUE)
  expect_match(
    sibling,
    "provider.x_row = np_conditional_cvls_provider_x_row;",
    fixed = TRUE
  )
  expect_match(sibling, "&provider", fixed = TRUE)
  expect_false(grepl("kernel_weighted_sum_np(", sibling, fixed = TRUE))
  expect_false(grepl("malloc(", sibling, fixed = TRUE))
  expect_false(grepl("calloc(", sibling, fixed = TRUE))
  expect_false(grepl("realloc(", sibling, fixed = TRUE))
})

test_that("conditional distribution route sibling delegates null and fails closed", {
  root <- locate_beta_conditional_cvls_sources()
  skip_if(is.null(root), "package sources unavailable")
  engine <- paste(
    readLines(file.path(root, "src", "jksum.c"), warn = FALSE),
    collapse = "\n"
  )
  sibling_start <- regexpr(
    "int np_conditional_distribution_cvls_lp_stream_ctx(",
    engine,
    fixed = TRUE
  )[[1L]]
  expect_gt(sibling_start, 0L)
  sibling <- substr(engine, sibling_start, nchar(engine))

  expect_match(
    sibling,
    paste0(
      "if(execution_context == NULL)\n",
      "    return np_conditional_distribution_cvls_lp_stream("
    ),
    fixed = TRUE
  )
  expect_match(
    sibling,
    "conditional distribution CVLS kernel route has an invalid layout",
    fixed = TRUE
  )
  expect_match(
    sibling,
    "provider.y_integral_row =\n    np_conditional_cvls_provider_y_integral_row;",
    fixed = TRUE
  )
  expect_match(
    sibling,
    "np_conditional_distribution_cvls_provider_supertile(",
    fixed = TRUE
  )
  expect_match(
    sibling,
    "np_conditional_distribution_cvls_provider_row_stream(",
    fixed = TRUE
  )
  expect_false(grepl("np_beta_objective_conditional_distribution_ls(",
                     sibling, fixed = TRUE))
})
