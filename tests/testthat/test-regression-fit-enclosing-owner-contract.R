locate_regression_fit_owner_source <- function() {
  roots <- unique(c(
    testthat::test_path("..", ".."),
    testthat::test_path("..", "..", ".."),
    Sys.getenv("R_PACKAGE_DIR", ""),
    Sys.getenv("R_PACKAGE_SOURCE", ""),
    getwd(),
    file.path(getwd(), "..")
  ))
  roots <- roots[nzchar(roots)]
  paths <- file.path(roots, "src", "jksum.c")
  paths <- paths[file.exists(paths)]
  if (!length(paths))
    return(NULL)
  paths[[1L]]
}

regression_fit_body <- function(source) {
  start <- grep(
    "^int kernel_estimate_regression_categorical_tree_np\\($", source
  )
  finish <- grep("^finish_regression_estimation:$", source)
  stopifnot(length(start) == 1L, length(finish) == 1L, finish > start)
  list(
    before_finish = paste(source[seq.int(start, finish - 1L)],
                          collapse = "\n"),
    after_finish = paste(source[seq.int(finish, finish + 70L)],
                         collapse = "\n")
  )
}

test_that("generic regression fit has one complete enclosing owner", {
  path <- locate_regression_fit_owner_source()
  skip_if(is.null(path), "package sources unavailable")
  source <- readLines(path, warn = FALSE)
  all_source <- paste(source, collapse = "\n")
  compact <- gsub("[[:space:]]+", " ", all_source)
  body <- regression_fit_body(source)
  before <- gsub("[[:space:]]+", " ", body$before_finish)
  after <- gsub("[[:space:]]+", " ", body$after_finish)

  for (resource in c(
    "operator", "kernel_c", "kernel_u", "kernel_o", "lambda",
    "matrix_bandwidth", "matrix_bandwidth_deriv", "ordered_tables",
    "matrix_ordered_indices", "ordered_indices_storage", "gate_context",
    "ov_cont_ok", "ov_cont_hmin", "ov_cont_k0", "ov_disc_uno_ok",
    "ov_disc_uno_const", "ov_disc_ord_ok", "ov_disc_ord_const"
  ))
    expect_match(compact, resource, fixed = TRUE, info = resource)

  expect_match(
    compact,
    paste0(
      "np_disc_profile_cache_clear(); np_cont_largeh_cache_clear(); ",
      "for(coordinate = 0; "
    ),
    fixed = TRUE
  )
  expect_match(
    compact,
    paste0(
      "free(owner->operator); free(owner->kernel_c); ",
      "free(owner->kernel_u); free(owner->kernel_o); free(owner->lambda);"
    ),
    fixed = TRUE
  )
  expect_match(
    compact,
    paste0(
      "free_tmat(owner->matrix_bandwidth); ",
      "free_mat(owner->matrix_bandwidth_deriv, owner->num_reg_continuous);"
    ),
    fixed = TRUE
  )
  expect_match(
    compact,
    paste0(
      "R_UnwindProtect( np_regression_fit_bandwidth_execute, call, ",
      "np_regression_fit_owner_cleanup, owner, NULL);"
    ),
    fixed = TRUE
  )
  expect_match(
    compact,
    "const NP_GateOverrideCtx *gate_context;",
    fixed = TRUE
  )
  expect_equal(
    sum(gregexpr(".enclosing_owner = &fit_owner", before,
                 fixed = TRUE)[[1L]] > 0L),
    3L
  )
  expect_match(
    compact,
    paste0(
      "if(jump && owner->enclosing_owner != NULL) ",
      "np_regression_fit_owner_clear(owner->enclosing_owner);"
    ),
    fixed = TRUE
  )
  expect_match(
    compact,
    paste0(
      "if(jump && call->enclosing_owner != NULL) ",
      "np_regression_fit_owner_clear(call->enclosing_owner);"
    ),
    fixed = TRUE
  )
  expect_match(
    compact,
    paste0(
      "np_progress_fit_loop_step_owned( i + 1, fit_progress_total, ",
      "np_regression_fit_owner_cleanup, (void *)&fit_owner);"
    ),
    fixed = TRUE
  )

  beta_return <- regexpr(
    "np_beta_scalar_regression_fit_canonical(", before, fixed = TRUE
  )[[1L]]
  refresh <- regexpr("np_refresh_runtime_tolerances();", before, fixed = TRUE)[[1L]]
  owner_init <- regexpr(
    "np_regression_fit_owner_init( &fit_owner", before, fixed = TRUE
  )[[1L]]
  bandwidth <- regexpr(
    "np_regression_fit_bandwidth_prepare_owned(", before, fixed = TRUE
  )[[1L]]
  expect_gt(beta_return, 0L)
  expect_gt(refresh, beta_return)
  expect_gt(owner_init, refresh)
  expect_gt(bandwidth, owner_init)

  owned_region <- substr(before, owner_init, nchar(before))
  expect_false(grepl("error(", owned_region, fixed = TRUE))
  expect_false(grepl("alloc_vecd(", owned_region, fixed = TRUE))
  expect_false(grepl("alloc_tmatd(", owned_region, fixed = TRUE))
  expect_false(grepl("alloc_matd(", owned_region, fixed = TRUE))
  expect_false(grepl("np_jksum_malloc_array3_or_die(", owned_region,
                     fixed = TRUE))
  expect_false(grepl("num_obs_train\\s*\\*\\s*num_obs_train", owned_region))
  expect_false(grepl("num_obs_eval\\s*\\*\\s*num_obs_eval", owned_region))

  expect_match(after, "np_regression_fit_owner_clear(&fit_owner);",
               fixed = TRUE)
  cleanup <- regexpr(
    "np_regression_fit_owner_clear(&fit_owner);", after, fixed = TRUE
  )[[1L]]
  first_error <- regexpr("error(", after, fixed = TRUE)[[1L]]
  expect_gt(cleanup, 0L)
  expect_gt(first_error, cleanup)
})

test_that("MPI owner-row solve exhaustion is not a beta-row failure", {
  path <- locate_regression_fit_owner_source()
  skip_if(is.null(path), "package sources unavailable")
  source <- paste(readLines(path, warn = FALSE), collapse = "\n")

  expect_match(
    source,
    "NP_REGRESSION_GENERAL_LP_FIT_ERR_OWNER_SOLVE = -9",
    fixed = TRUE
  )
  expect_match(
    source,
    'error("LP solve failed in MPI owner-row path")',
    fixed = TRUE
  )
  if (grepl("#define MPI2", source, fixed = TRUE) ||
      grepl("NPRegMpiOwnerChunk mpi_owner_chunk;", source, fixed = TRUE)) {
    expect_match(
      source,
      "execution->status = NP_REGRESSION_GENERAL_LP_FIT_ERR_OWNER_SOLVE;",
      fixed = TRUE
    )
  }
})

