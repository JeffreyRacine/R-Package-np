np_conditional_density_prepared_source <- function() {
  source_file <- npRmpi_test_source_path("src", "np.c")
  paste(readLines(source_file, warn = FALSE), collapse = "\n")
}

np_conditional_density_prepared_function <- function(source, name,
                                                      return_type = "void") {
  markers <- gregexpr(
    sprintf("(?:static\\s+)?%s\\s+%s\\s*\\(", return_type, name),
    source, perl = TRUE
  )[[1L]]
  expect_false(identical(markers, -1L))
  for (marker in markers) {
    tail <- substring(source, marker)
    open <- regexpr("\\{", tail, perl = TRUE)[[1L]]
    semi <- regexpr(";", tail, fixed = TRUE)[[1L]]
    if (open <= 0L || (semi > 0L && semi < open))
      next
    chars <- strsplit(substring(tail, open), "", fixed = TRUE)[[1L]]
    depth <- 0L
    for (i in seq_along(chars)) {
      if (chars[[i]] == "{")
        depth <- depth + 1L
      else if (chars[[i]] == "}") {
        depth <- depth - 1L
        if (depth == 0L)
          return(substring(tail, 1L, open + i - 1L))
      }
    }
  }
  fail(sprintf("could not extract %s", name))
}

test_that("MPI conditional-density bandwidth state has one typed owner", {
  source <- np_conditional_density_prepared_source()
  clear <- np_conditional_density_prepared_function(
    source, "np_conditional_density_prepared_context_clear_internal"
  )

  expect_match(source, "} NPConditionalDensityPreparedCtx;", fixed = TRUE)
  expect_equal(
    lengths(regmatches(
      source,
      gregexpr(
        "static NPConditionalDensityPreparedCtx np_conditional_density_prepared_context =",
        source, fixed = TRUE
      )
    )),
    1L
  )
  for (field in c(
    "matrix_y_unordered_original", "matrix_y_ordered_original",
    "matrix_y_continuous_original", "matrix_xy_unordered_original",
    "matrix_xy_ordered_original", "matrix_xy_continuous_original",
    "extendednn_upper", "powell_directions", "vector_scale_factor",
    "scale_factor_startbest", "powell_step", "ipt_x", "ipt_y", "ipt_xy",
    "penalty_ready", "penalty_value"
  )) {
    expect_match(
      clear,
      paste0("np_conditional_density_prepared_context.", field),
      fixed = TRUE, info = field
    )
  }
  expect_match(
    source, "np_conditional_density_prepared_context_prepare_internal(",
    fixed = TRUE
  )
  expect_match(
    source, "np_conditional_density_prepared_context_clear_internal();",
    fixed = TRUE
  )
  expect_match(clear, "bwm_clear_floor_context();", fixed = TRUE)
  expect_match(
    source, "np_conditional_density_refresh_penalty_canonical(", fixed = TRUE
  )
  expect_match(
    source,
    "bwm_penalty_mode = np_conditional_density_prepared_context.penalty_ready;",
    fixed = TRUE
  )
  expect_false(grepl("nomad_shadow", source, fixed = TRUE))
})

test_that("MPI conditional-density searches share bounded prepared scratch", {
  source <- np_conditional_density_prepared_source()
  evaluator <- np_conditional_density_prepared_function(
    source, "np_conditional_density_prepared_context_eval_native_raw", "int"
  )

  expect_match(evaluator, ".eval_bandwidth", fixed = TRUE)
  expect_match(evaluator, ".eval_degree", fixed = TRUE)
  expect_false(grepl("NP_NOMAD_CALLBACK_CALLOC", evaluator, fixed = TRUE))
  expect_false(grepl("NP_NOMAD_CALLBACK_FREE", evaluator, fixed = TRUE))
  expect_false(grepl(
    "goto canonical_conditional_density_search;", source, fixed = TRUE
  ))
})
