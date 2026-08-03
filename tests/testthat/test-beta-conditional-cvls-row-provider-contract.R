locate_beta_conditional_cvls_provider_source <- function() {
  roots <- unique(c(
    Sys.getenv("NP_SOURCE_ROOT", unset = ""),
    normalizePath(file.path(getwd(), "..", ".."), mustWork = FALSE),
    normalizePath(getwd(), mustWork = FALSE)
  ))
  roots <- roots[nzchar(roots)]
  for (root in roots) {
    path <- file.path(root, "src", "jksum.c")
    if (file.exists(path))
      return(path)
  }
  NULL
}

test_that("bounded conditional CVLS has one typed row-provider seam", {
  path <- locate_beta_conditional_cvls_provider_source()
  skip_if(is.null(path), "package source unavailable")
  engine <- paste(readLines(path, warn = FALSE), collapse = "\n")

  expect_match(engine, "} NPConditionalCVLSRowProvider;", fixed = TRUE)
  for (member in c("x_row", "y_train_row", "y_convolution_row",
                   "y_eval_block", "y_integral_row"))
    expect_match(engine, paste0("(*", member, ")"), fixed = TRUE)
  expect_match(
    engine,
    paste0(
      "np_conditional_density_cvls_bounded_i1_quadrature_row_stream(",
      "vector_scale_factor,\n",
      "                                                                        cv,\n",
      "                                                                        NULL);"
    ),
    fixed = TRUE
  )
  expect_match(
    engine,
    paste0(
      "np_conditional_density_cvls_bounded_i1_quadrature_general_row_stream(",
      "vector_scale_factor,\n",
      "                                                                                cv,\n",
      "                                                                                NULL);"
    ),
    fixed = TRUE
  )
  expect_false(grepl("NP_BOUNDED_CVLS_I1_MODE_", engine, fixed = TRUE))
  expect_false(grepl("i1_mode", engine, fixed = TRUE))
})

test_that("conditional distribution provider is bounded and operator explicit", {
  path <- locate_beta_conditional_cvls_provider_source()
  skip_if(is.null(path), "package source unavailable")
  engine <- paste(readLines(path, warn = FALSE), collapse = "\n")
  start <- regexpr(
    "static int np_conditional_distribution_cvls_provider_supertile(",
    engine,
    fixed = TRUE
  )[[1L]]
  finish <- regexpr(
    "static int np_conditional_distribution_cvls_lp_row_stream(",
    engine,
    fixed = TRUE
  )[[1L]]
  expect_gt(start, 0L)
  expect_gt(finish, start)
  route <- substr(engine, start, finish - 1L)

  expect_match(route, "group_width = MIN(4, MAX(1, nblocks));",
               fixed = TRUE)
  expect_match(route, "np_blas_dgemm_tn_int(", fixed = TRUE)
  expect_match(route, "provider->y_integral_row(", fixed = TRUE)
  expect_match(route, "np_continuous_kernel_scaled_restore(", fixed = TRUE)
  expect_false(grepl(
    "np_cvls_workspace_matrix_try(num_train, num_train",
    route,
    fixed = TRUE
  ))
  expect_false(grepl("np_beta_objective_conditional_distribution_ls(",
                     route, fixed = TRUE))
  expect_match(
    engine,
    paste0(
      "&route_context, vector_scale_factor, execution_context,\n",
      "       OP_INTEGRAL)"
    ),
    fixed = TRUE
  )
})

test_that("live bounded CVLS retains literal legacy arithmetic under null", {
  path <- locate_beta_conditional_cvls_provider_source()
  skip_if(is.null(path), "package source unavailable")
  engine <- paste(readLines(path, warn = FALSE), collapse = "\n")

  expect_match(
    engine,
    "if(provider == NULL) {\n    if(np_conditional_xrow_ctx_prepare(",
    fixed = TRUE
  )
  expect_match(
    engine,
    "if(np_conditional_y_scalar_eval_from_ctx(",
    fixed = TRUE
  )
  expect_match(
    engine,
    "if(provider->y_eval_block(provider->context, q,",
    fixed = TRUE
  )
  expect_match(
    engine,
    "} else if(np_conditional_y_eval_any_block_stream_core(",
    fixed = TRUE
  )
  expect_match(
    engine,
    "np_conditional_x_weight_row_stream_core(",
    fixed = TRUE
  )
  expect_match(engine, "np_conditional_yrow_from_ctx(", fixed = TRUE)
  expect_match(
    engine,
    paste0(
      "np_conditional_density_cvls_lp_row_stream(\n",
      "      vector_scale_factor, cv, NULL);"
    ),
    fixed = TRUE
  )

  sibling_start <- regexpr(
    "int np_conditional_density_cvls_lp_stream_ctx(", engine, fixed = TRUE
  )[[1L]]
  expect_gt(sibling_start, 0L)
  sibling <- substr(engine, sibling_start, nchar(engine))
  expect_match(sibling, "NPConditionalCVLSRowProvider provider;", fixed = TRUE)
  expect_match(sibling, "provider.y_eval_block =", fixed = TRUE)
  expect_match(sibling, "&provider", fixed = TRUE)
})

test_that("analytic routed CVLS uses bounded tiles and only an allocation fallback", {
  path <- locate_beta_conditional_cvls_provider_source()
  skip_if(is.null(path), "package source unavailable")
  engine <- paste(readLines(path, warn = FALSE), collapse = "\n")
  start <- regexpr(
    "static int np_conditional_density_cvls_provider_supertile_stream(",
    engine,
    fixed = TRUE
  )[[1L]]
  finish <- regexpr(
    "int np_conditional_density_cvls_lp_stream_ctx(",
    engine,
    fixed = TRUE
  )[[1L]]
  expect_gt(start, 0L)
  expect_gt(finish, start)
  route <- substr(engine, start, finish - 1L)

  expect_match(
    route,
    "np_conditional_lp_cvls_block_size(num_obs, 6U, 0U)",
    fixed = TRUE
  )
  expect_match(route, "group_width = MIN(4, MAX(1, nblocks));", fixed = TRUE)
  expect_match(route, "np_blas_dgemm_tn_int(", fixed = TRUE)
  expect_false(grepl(
    "np_cvls_workspace_matrix_try(num_obs, num_obs",
    route,
    fixed = TRUE
  ))
  expect_false(grepl("kernel_weighted_sum_np(", route, fixed = TRUE))

  sibling <- substr(engine, finish, nchar(engine))
  expect_match(
    sibling,
    paste0(
      "if(status == 2)\n",
      "      status = np_conditional_density_cvls_lp_row_stream("
    ),
    fixed = TRUE
  )
})
