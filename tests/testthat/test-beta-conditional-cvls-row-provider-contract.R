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

extract_beta_conditional_cvls_region <- function(engine,
                                                 start_marker,
                                                 finish_marker = NULL,
                                                 start_occurrence = c("first", "last")) {
  start_occurrence <- match.arg(start_occurrence)
  starts <- gregexpr(start_marker, engine, fixed = TRUE)[[1L]]
  starts <- starts[starts > 0L]
  start <- if (identical(start_occurrence, "last")) {
    tail(starts, 1L)
  } else {
    starts[[1L]]
  }
  finish <- if (is.null(finish_marker)) {
    nchar(engine) + 1L
  } else {
    finishes <- gregexpr(finish_marker, engine, fixed = TRUE)[[1L]]
    finishes <- finishes[finishes > start]
    finishes[[1L]]
  }

  expect_gt(start, 0L)
  expect_gt(finish, start)
  substr(engine, start, finish - 1L)
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

test_that("local conditional distribution supertile remains ownership-free", {
  path <- locate_beta_conditional_cvls_provider_source()
  skip_if(is.null(path), "package source unavailable")
  engine <- paste(readLines(path, warn = FALSE), collapse = "\n")
  route <- extract_beta_conditional_cvls_region(
    engine,
    "static int np_conditional_distribution_cvls_provider_supertile(",
    "static int np_conditional_distribution_cvls_lp_row_stream("
  )

  expect_match(route, "group_width = MIN(4, MAX(1, nblocks));",
               fixed = TRUE)
  expect_match(route, "np_blas_dgemm_tn_int(", fixed = TRUE)
  expect_match(route, "provider->y_integral_row(", fixed = TRUE)
  expect_match(route, "np_continuous_kernel_scaled_restore(", fixed = TRUE)
  expect_false(grepl("np_objective_outer_", route, fixed = TRUE))
  expect_false(grepl("contributions", route, fixed = TRUE))
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

test_that("local conditional distribution row stream remains ownership-free", {
  path <- locate_beta_conditional_cvls_provider_source()
  skip_if(is.null(path), "package source unavailable")
  engine <- paste(readLines(path, warn = FALSE), collapse = "\n")
  route <- extract_beta_conditional_cvls_region(
    engine,
    "static int np_conditional_distribution_cvls_provider_row_stream(",
    "static int np_conditional_distribution_cvls_provider_supertile("
  )

  expect_match(route, "provider->x_row(", fixed = TRUE)
  expect_match(route, "provider->y_integral_row(", fixed = TRUE)
  expect_match(route, "np_blas_ddot_int(", fixed = TRUE)
  expect_false(grepl("np_objective_outer_", route, fixed = TRUE))
  expect_false(grepl("contributions", route, fixed = TRUE))
  expect_false(grepl("alloc_matd(num_train, num_train)", route, fixed = TRUE))
  expect_false(grepl("diag(num_train)", route, fixed = TRUE))
})

test_that("parallel conditional distribution row fallback retains ownership", {
  path <- locate_beta_conditional_cvls_provider_source()
  skip_if(is.null(path), "package source unavailable")
  engine <- paste(readLines(path, warn = FALSE), collapse = "\n")
  route <- extract_beta_conditional_cvls_region(
    engine,
    "static int np_conditional_distribution_cvls_provider_row_stream_parallel(",
    "static int np_conditional_distribution_cvls_provider_supertile_parallel(",
    start_occurrence = "last"
  )

  expect_match(
    route,
    paste0(
      "np_objective_outer_rows_enabled(\n",
      "       int_conditional_prepared_context_extern)"
    ),
    fixed = TRUE
  )
  expect_match(route, "np_objective_outer_buffer_prepare(", fixed = TRUE)
  expect_match(route, "np_objective_outer_owned_rows(", fixed = TRUE)
  expect_match(route, "contributions[i] = row_loss;", fixed = TRUE)
  expect_match(route, "np_objective_outer_buffer_finish(", fixed = TRUE)
  expect_match(route, '"NP_RMPI_INJECT_CDIST_CVLS_FAIL_RANK"', fixed = TRUE)
  expect_false(grepl("alloc_matd(num_train, num_train)", route, fixed = TRUE))
  expect_false(grepl("diag(num_train)", route, fixed = TRUE))
})

test_that("parallel conditional distribution supertile owns bounded groups", {
  path <- locate_beta_conditional_cvls_provider_source()
  skip_if(is.null(path), "package source unavailable")
  engine <- paste(readLines(path, warn = FALSE), collapse = "\n")
  route <- extract_beta_conditional_cvls_region(
    engine,
    "static int np_conditional_distribution_cvls_provider_supertile_parallel(",
    start_occurrence = "last"
  )

  expect_match(
    route,
    "nblocks/owned_stride + ((nblocks % owned_stride) != 0)",
    fixed = TRUE
  )
  expect_match(route, "np_objective_outer_buffer_prepare(", fixed = TRUE)
  expect_match(route, "const int first_owned_block = my_rank;", fixed = TRUE)
  expect_match(route, "block_index[g] = i0 + g*owned_stride;", fixed = TRUE)
  expect_match(route, "contributions[block_index[g]] = block_sum[g];",
               fixed = TRUE)
  expect_match(route, "np_objective_outer_buffer_finish(", fixed = TRUE)
  expect_match(route, '"NP_RMPI_INJECT_CDIST_CVLS_FAIL_RANK"', fixed = TRUE)
  expect_match(route, "np_blas_dgemm_tn_int(", fixed = TRUE)
  expect_match(route, "provider->y_integral_row(", fixed = TRUE)
  expect_match(route, "np_continuous_kernel_scaled_restore(", fixed = TRUE)
  expect_false(grepl(
    "np_cvls_workspace_matrix_try(num_train, num_train",
    route,
    fixed = TRUE
  ))
  expect_false(grepl("alloc_matd(num_train, num_train)", route, fixed = TRUE))
  expect_false(grepl("diag(num_train)", route, fixed = TRUE))
})

test_that("conditional distribution dispatcher isolates local and MPI owners", {
  path <- locate_beta_conditional_cvls_provider_source()
  skip_if(is.null(path), "package source unavailable")
  engine <- paste(readLines(path, warn = FALSE), collapse = "\n")
  route <- substr(
    engine,
    regexpr(
      "int np_conditional_distribution_cvls_lp_stream_ctx(",
      engine,
      fixed = TRUE
    )[[1L]],
    nchar(engine)
  )

  expect_match(
    route,
    paste0(
      "if(np_objective_outer_rows_enabled(\n",
      "       int_conditional_prepared_context_extern)) {"
    ),
    fixed = TRUE
  )
  expect_match(
    route,
    "np_conditional_distribution_cvls_provider_supertile_parallel(",
    fixed = TRUE
  )
  expect_match(
    route,
    "np_conditional_distribution_cvls_provider_row_stream_parallel(",
    fixed = TRUE
  )
  expect_match(
    route,
    "np_conditional_distribution_cvls_provider_supertile(",
    fixed = TRUE
  )
  expect_match(
    route,
    "np_conditional_distribution_cvls_provider_row_stream(",
    fixed = TRUE
  )
  expect_false(grepl("np_beta_objective_conditional_distribution_ls(",
                     route, fixed = TRUE))
})

test_that("conditional distribution preparation arms the shared owner flag", {
  path <- locate_beta_conditional_cvls_provider_source()
  skip_if(is.null(path), "package source unavailable")
  native_path <- file.path(dirname(path), "np.c")
  skip_if_not(file.exists(native_path), "package np.c source unavailable")
  native <- paste(readLines(native_path, warn = FALSE), collapse = "\n")

  starts <- gregexpr(
    "static void np_conditional_distribution_prepared_context_destroy(",
    native,
    fixed = TRUE
  )[[1L]]
  prepares <- gregexpr(
    "static void np_distribution_conditional_bw_mode(",
    native,
    fixed = TRUE
  )[[1L]]
  start <- tail(starts[starts > 0L], 1L)
  prepare <- tail(prepares[prepares > 0L], 1L)
  expect_gt(start, 0L)
  expect_gt(prepare, start)
  owner <- substr(native, start, prepare - 1L)
  expect_match(
    owner,
    "int_conditional_prepared_context_extern = 0;",
    fixed = TRUE
  )

  expect_match(
    native,
    paste0(
      "prepared_context->active = 1;\n",
      "  int_conditional_prepared_context_extern = 1;\n",
      "  prepared_context->cdfontrain = cdfontrain;"
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
