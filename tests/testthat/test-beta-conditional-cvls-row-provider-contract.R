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
                   "y_scalar_eval_row", "y_eval_block"))
    expect_match(engine, paste0("(*", member, ")"), fixed = TRUE)
  expect_match(
    engine,
    paste0(
      "np_conditional_density_cvls_bounded_i1_quadrature_row_stream(",
      "vector_scale_factor,\n",
      "                                                                        cv,\n",
      "                                                                        NP_BOUNDED_CVLS_I1_MODE_BOOK,\n",
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
      "                                                                                NP_BOUNDED_CVLS_I1_MODE_BOOK,\n",
      "                                                                                NULL);"
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
    "} else if(np_conditional_y_scalar_eval_from_ctx(",
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
  expect_match(sibling, "return 1;", fixed = TRUE)
  expect_false(grepl("NPConditionalCVLSRowProvider provider", sibling,
                     fixed = TRUE))
})
