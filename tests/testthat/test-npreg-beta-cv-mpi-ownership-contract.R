npreg_beta_cv_mpi_source <- function() {
  candidates <- c(
    test_path("..", "..", "src", "jksum.c"),
    test_path("..", "..", "..", "src", "jksum.c"),
    file.path(Sys.getenv("R_PACKAGE_DIR", ""), "src", "jksum.c"),
    file.path(Sys.getenv("R_PACKAGE_SOURCE", ""), "src", "jksum.c"),
    file.path(getwd(), "src", "jksum.c"),
    file.path(getwd(), "..", "src", "jksum.c")
  )
  candidates <- unique(candidates[nzchar(candidates)])
  hits <- candidates[file.exists(candidates)]
  if (!length(hits)) return(NULL)
  paste(readLines(hits[[1L]], warn = FALSE), collapse = "\n")
}

npreg_beta_cv_mpi_body <- function(source, start, stop) {
  first <- regexpr(start, source, fixed = TRUE)[[1L]]
  expect_gt(first, 0L)
  remainder <- substr(source, first, nchar(source))
  offset <- regexpr(stop, remainder, fixed = TRUE)[[1L]]
  expect_gt(offset, 0L)
  substr(source, first, first + offset - 2L)
}

test_that("fixed scalar CVLS admits one complete-row MPI owner", {
  source <- npreg_beta_cv_mpi_source()
  skip_if(is.null(source), "package C source unavailable")

  body <- npreg_beta_cv_mpi_body(
    source,
    "static int np_regression_cv_scalar_continuous_route_parallel_body(",
    "static SEXP np_regression_cv_scalar_route_execute"
  )
  expect_match(
    body,
    "call->bwm != RBWM_CVLS",
    fixed = TRUE
  )
  expect_match(
    body,
    "call->bandwidth_mode != BW_FIXED",
    fixed = TRUE
  )
  expect_match(body, "np_objective_outer_buffer_prepare(", fixed = TRUE)
  expect_match(body, "np_objective_outer_owned_rows(", fixed = TRUE)
  expect_match(body, "contributions[evaluation] = row_loss;", fixed = TRUE)
  expect_match(body, "np_objective_outer_buffer_finish(", fixed = TRUE)
  expect_match(
    body,
    '"NP_RMPI_INJECT_REG_ROUTED_CV_FAIL_RANK"',
    fixed = TRUE
  )
  expect_match(
    body,
    "for(evaluation = 0; evaluation < num_obs; ++evaluation) {\n    cv += contributions[evaluation];",
    fixed = TRUE
  )
  expect_false(grepl("evaluation % iNum_Processors", body, fixed = TRUE))
  expect_false(grepl("MPI_Allreduce(", body, fixed = TRUE))
  expect_false(grepl("alloc_matd(num_obs, num_obs)", body, fixed = TRUE))
  expect_false(grepl("diag(num_obs)", body, fixed = TRUE))
})

test_that("scalar owner admits both fixed objectives and retains NN isolation", {
  source <- npreg_beta_cv_mpi_source()
  skip_if(is.null(source), "package C source unavailable")

  body <- npreg_beta_cv_mpi_body(
    source,
    "static SEXP np_regression_cv_scalar_route_execute",
    "static NP_NOINLINE int np_regression_cv_scalar_continuous_route("
  )
  gate <- regmatches(
    body,
    regexpr(
      "np_objective_outer_rows_enabled\\([^;]+;",
      body,
      perl = TRUE
    )
  )
  expect_length(gate, 1L)
  expect_match(gate, "RBWM_CVLS", fixed = TRUE)
  expect_match(gate, "RBWM_CVAIC", fixed = TRUE)
  expect_match(gate, "BW_FIXED", fixed = TRUE)
  expect_false(grepl("BW_GEN_NN", gate, fixed = TRUE))
  expect_false(grepl("BW_ADAP_NN", gate, fixed = TRUE))
})

test_that("fixed wider LP CVLS admits one complete-row MPI owner", {
  source <- npreg_beta_cv_mpi_source()
  skip_if(is.null(source), "package C source unavailable")

  body <- npreg_beta_cv_mpi_body(
    source,
    "static int np_regression_cv_lp_continuous_route_parallel_body(",
    "static SEXP np_regression_cv_lp_route_execute"
  )
  expect_match(body, "call->bwm != RBWM_CVLS", fixed = TRUE)
  expect_match(body, "call->bandwidth_mode != BW_FIXED", fixed = TRUE)
  expect_match(body, "np_objective_outer_buffer_prepare(", fixed = TRUE)
  expect_match(body, "np_objective_outer_owned_rows(", fixed = TRUE)
  expect_match(
    body,
    "contributions[evaluation] = residual*residual;",
    fixed = TRUE
  )
  expect_match(body, "np_objective_outer_buffer_finish(", fixed = TRUE)
  expect_match(
    body,
    '"NP_RMPI_INJECT_REG_ROUTED_CV_FAIL_RANK"',
    fixed = TRUE
  )
  expect_match(
    body,
    "for(evaluation = 0; evaluation < num_obs; ++evaluation) {\n    cv += contributions[evaluation];",
    fixed = TRUE
  )
  expect_false(grepl("RBWM_CVAIC", body, fixed = TRUE))
  expect_false(grepl("evaluation % iNum_Processors", body, fixed = TRUE))
  expect_false(grepl("MPI_Allreduce(", body, fixed = TRUE))
  expect_false(grepl("alloc_matd(num_obs, num_obs)", body, fixed = TRUE))
  expect_false(grepl("diag(num_obs)", body, fixed = TRUE))
})

test_that("wider LP owner admits fixed CVLS only", {
  source <- npreg_beta_cv_mpi_source()
  skip_if(is.null(source), "package C source unavailable")

  body <- npreg_beta_cv_mpi_body(
    source,
    "static SEXP np_regression_cv_lp_route_execute",
    "static NP_NOINLINE int np_regression_cv_lp_continuous_route("
  )
  gate <- regmatches(
    body,
    regexpr(
      "np_objective_outer_rows_enabled\\([^;]+;",
      body,
      perl = TRUE
    )
  )
  expect_length(gate, 1L)
  expect_match(gate, "RBWM_CVLS", fixed = TRUE)
  expect_match(gate, "BW_FIXED", fixed = TRUE)
  expect_false(grepl("RBWM_CVAIC", gate, fixed = TRUE))
  expect_false(grepl("BW_GEN_NN", gate, fixed = TRUE))
  expect_false(grepl("BW_ADAP_NN", gate, fixed = TRUE))
})

test_that("incumbent local wider LP hot loop remains ownership-free", {
  source <- npreg_beta_cv_mpi_source()
  skip_if(is.null(source), "package C source unavailable")

  body <- npreg_beta_cv_mpi_body(
    source,
    "static int np_regression_cv_lp_continuous_route_body(",
    "static int np_regression_cv_lp_continuous_route_parallel_body("
  )
  expect_false(grepl("np_objective_outer_", body, fixed = TRUE))
  expect_false(grepl("contributions", body, fixed = TRUE))
  expect_match(
    body,
    "for(evaluation = 0; evaluation < num_obs; ++evaluation)",
    fixed = TRUE
  )
})

test_that("incumbent local scalar hot loop remains ownership-free", {
  source <- npreg_beta_cv_mpi_source()
  skip_if(is.null(source), "package C source unavailable")

  body <- npreg_beta_cv_mpi_body(
    source,
    "static int np_regression_cv_scalar_continuous_route_body(",
    "static int np_regression_cv_scalar_continuous_route_parallel_body("
  )
  expect_false(grepl("np_objective_outer_", body, fixed = TRUE))
  expect_false(grepl("contributions", body, fixed = TRUE))
  expect_match(
    body,
    "for(evaluation = 0; evaluation < num_obs; ++evaluation)",
    fixed = TRUE
  )
  expect_match(
    body,
    "&cv, &traceH) != 0",
    fixed = TRUE
  )
})

test_that("parallel rows reuse the incumbent scalar estimator finisher", {
  source <- npreg_beta_cv_mpi_source()
  skip_if(is.null(source), "package C source unavailable")

  helper <- npreg_beta_cv_mpi_body(
    source,
    "static inline int np_regression_cv_scalar_accumulate_scaled_row(",
    "typedef struct {\n  int bwm;"
  )
  expect_match(
    helper,
    "numerator += row[observation]*response[observation];",
    fixed = TRUE
  )
  expect_match(helper, "fitted = numerator/row_sum;", fixed = TRUE)
  expect_match(helper, "residual = response[evaluation] - fitted;",
               fixed = TRUE)
  expect_match(helper, "*cv += residual*residual;", fixed = TRUE)
  expect_match(
    helper,
    "*traceH += row[evaluation]/row_sum;",
    fixed = TRUE
  )
  expect_match(
    source,
    "&row_loss, &row_trace) != 0",
    fixed = TRUE
  )
  expect_match(
    source,
    "contributions[num_obs + evaluation] = row_trace;",
    fixed = TRUE
  )
})

test_that("local and parallel scalar routes share the objective finisher", {
  source <- npreg_beta_cv_mpi_source()
  skip_if(is.null(source), "package C source unavailable")

  helper <- npreg_beta_cv_mpi_body(
    source,
    "static NP_ALWAYS_INLINE int np_regression_cv_finish_objective(",
    "typedef struct {\n  int bwm;"
  )
  expect_match(helper, "cv /= (double)num_obs;", fixed = TRUE)
  expect_match(helper, "if(bwm == RBWM_CVAIC)", fixed = TRUE)
  expect_match(helper, "cv = log(cv) + penalty;", fixed = TRUE)
  expect_match(
    source,
    "bwm, num_obs, cv, traceH, objective",
    fixed = TRUE
  )
  expect_match(
    source,
    "call->bwm, num_obs, cv, traceH, call->objective",
    fixed = TRUE
  )
  expect_match(
    source,
    "bwm, num_obs, cv, traceH, objective",
    fixed = TRUE
  )
  expect_match(
    source,
    "RBWM_CVLS, num_obs, cv, 0.0, call->objective",
    fixed = TRUE
  )
})
