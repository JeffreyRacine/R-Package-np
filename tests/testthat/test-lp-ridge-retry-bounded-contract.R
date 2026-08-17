locate_lp_retry_source <- function(file) {
  candidates <- c(
    test_path("..", "..", "src", file),
    test_path("..", "..", "..", "src", file),
    file.path(Sys.getenv("R_PACKAGE_DIR", ""), "src", file),
    file.path(Sys.getenv("R_PACKAGE_SOURCE", ""), "src", file),
    file.path(getwd(), "src", file),
    file.path(getwd(), "..", "src", file)
  )
  candidates <- unique(candidates[nzchar(candidates)])
  hits <- candidates[file.exists(candidates)]
  if (!length(hits))
    return(NULL)
  hits[[1L]]
}

locate_lp_retry_r_source <- function(file) {
  candidates <- c(
    test_path("..", "..", "R", file),
    test_path("..", "..", "..", "R", file),
    file.path(Sys.getenv("R_PACKAGE_DIR", ""), "R", file),
    file.path(Sys.getenv("R_PACKAGE_SOURCE", ""), "R", file),
    file.path(getwd(), "R", file),
    file.path(getwd(), "..", "R", file)
  )
  candidates <- unique(candidates[nzchar(candidates)])
  hits <- candidates[file.exists(candidates)]
  if (!length(hits))
    return(NULL)
  hits[[1L]]
}

test_that("all canonical LP solve retries are bounded", {
  jksum.file <- locate_lp_retry_source("jksum.c")
  solve.c.file <- locate_lp_retry_source("jksum_lp_solve.c")
  solve.h.file <- locate_lp_retry_source("jksum_lp_solve.h")
  reghat.c.file <- locate_lp_retry_source("reghat_fast.c")
  reghat.r.file <- locate_lp_retry_r_source("np.reghat.R")
  skip_if(any(vapply(
    list(jksum.file, solve.c.file, solve.h.file, reghat.c.file, reghat.r.file),
    is.null,
    logical(1L)
  )), "complete LP retry sources unavailable in this test context")

  jksum <- paste(readLines(jksum.file, warn = FALSE), collapse = "\n")
  solve.lines <- readLines(solve.c.file, warn = FALSE)
  solve.c <- paste(solve.lines, collapse = "\n")
  solve.h <- paste(readLines(solve.h.file, warn = FALSE), collapse = "\n")
  reghat.c <- paste(readLines(reghat.c.file, warn = FALSE), collapse = "\n")
  reghat.r <- paste(readLines(reghat.r.file, warn = FALSE), collapse = "\n")

  retry.count <- lengths(regmatches(
    jksum,
    gregexpr(
      "while\\s*\\(\\s*!np_lp_solve_workspace_solve\\(",
      jksum
    )
  ))
  finite.guard.count <- lengths(regmatches(
    jksum,
    gregexpr("np_lp_solve_workspace_sources_finite\\(", jksum)
  ))
  step.guard.count <- lengths(regmatches(
    jksum,
    gregexpr("ridge_steps >= NP_LP_SOLVE_MAX_RIDGE_STEPS", jksum)
  ))
  expect_identical(retry.count, 0L)
  expect_identical(finite.guard.count, 0L)
  expect_identical(step.guard.count, 0L)
  policy.body <- npRmpi_test_extract_c_function(
    solve.lines, "np_lp_solve_workspace_prepare_policy_factor"
  )
  expect_match(policy.body, "while(!np_lp_solve_workspace_factor(",
               fixed = TRUE)
  expect_match(policy.body, "np_lp_solve_workspace_sources_finite(",
               fixed = TRUE)
  expect_match(policy.body, "ridge_steps >= NP_LP_SOLVE_MAX_RIDGE_STEPS",
               fixed = TRUE)
  expect_true(grepl(
    "execution->status = NP_REGRESSION_GENERAL_LP_FIT_ERR_SOLVE;",
    jksum,
    fixed = TRUE
  ))
  expect_false(grepl("cleanup_glp_fit:", jksum, fixed = TRUE))
  expect_true(grepl(
    "#define NP_LP_SOLVE_MAX_RIDGE_STEPS 128",
    solve.h,
    fixed = TRUE
  ))
  expect_true(grepl(
    "np_lp_solve_workspace_sources_finite(",
    solve.c,
    fixed = TRUE
  ))
  expect_true(grepl(
    "int np_lp_solve_workspace_sources_finite(",
    solve.h,
    fixed = TRUE
  ))
  expect_true(grepl("NP_LP_COLD_PATH", solve.h, fixed = TRUE))
  expect_true(grepl(
    "for(i = 0; i < gram_elements; i++)",
    solve.c,
    fixed = TRUE
  ))
  expect_true(grepl(
    "for(i = 0; i < rhs_elements; i++)",
    solve.c,
    fixed = TRUE
  ))
  solve.body <- npRmpi_test_extract_c_function(
    solve.lines, "np_lp_solve_workspace_solve"
  )
  expect_false(grepl(
    "np_lp_solve_workspace_shape(",
    solve.body,
    fixed = TRUE
  ))
  expect_true(grepl(
    "(workspace == NULL) || (p <= 0) || (nrhs <= 0)",
    solve.body,
    fixed = TRUE
  ))

  expect_false(grepl("NP_LP_SOLVE_MAX_RIDGE_STEPS", reghat.c, fixed = TRUE))
  expect_false(grepl("np_reghat_sources_finite", reghat.c, fixed = TRUE))
  expect_true(grepl(
    "np_lp_solve_workspace_solve_adjoint(solve_workspace,",
    reghat.c,
    fixed = TRUE
  ))
  expect_false(grepl(
    "while\\s*\\(\\s*!np_reghat_solve_system\\(",
    reghat.c
  ))

  fallback.start <- regexpr(
    "\\.npreghat_exact_lp_matrix_from_kernel_weights <- function",
    reghat.r
  )
  fallback.stops <- c(
    regexpr(
      "\\.npreghat_exact_lp_apply_chunked_from_kernel_weights <- function",
      reghat.r
    ),
    regexpr(
      "\\.npreghat_exact_lp_apply_from_regression_core <- function",
      reghat.r
    )
  )
  fallback.stop <- min(fallback.stops[fallback.stops > fallback.start])
  expect_gt(fallback.start, 0L)
  expect_gt(fallback.stop, fallback.start)
  fallback <- substr(reghat.r, fallback.start, fallback.stop - 1L)
  expect_true(grepl("seq_len(128L)", fallback, fixed = TRUE))
  expect_true(grepl("non-finite system", fallback, fixed = TRUE))
  expect_false(grepl("repeat\\s*\\{", fallback))
})

test_that("compiled LP hat path rejects a non-finite system promptly", {
  package <- if ("npRmpi" %in% loadedNamespaces()) "npRmpi" else "np"
  kw <- matrix(1.0, nrow = 4L, ncol = 1L)
  wtrain <- cbind(1.0, c(NaN, 0.0, 0.0, 0.0))
  weval <- matrix(c(1.0, 0.0), nrow = 1L)

  expect_error(
    .Call(
      "C_np_reghat_lp_matrix_fast",
      kw,
      wtrain,
      weval,
      PACKAGE = package
    ),
    "LP solve failed in compiled hat-matrix path: non-finite system"
  )
})

test_that("MPI LP owner failures remain collective and fully owned", {
  jksum.file <- locate_lp_retry_source("jksum.c")
  skip_if(is.null(jksum.file), "source file src/jksum.c unavailable")
  source <- paste(readLines(jksum.file, warn = FALSE), collapse = "\n")

  expect_true(grepl(
    "owner->mpi_owner_chunk.sendbuf[0] = NAN;",
    source,
    fixed = TRUE
  ))
  expect_true(grepl(
    "np_reg_mpi_owner_chunk_allgather(&owner->mpi_owner_chunk);",
    source,
    fixed = TRUE
  ))
  expect_true(grepl(
    "!isfinite(np_reg_mpi_owner_chunk_recv_ptr(",
    source,
    fixed = TRUE
  ))
  free.pos <- regexpr(
    "np_reg_mpi_owner_chunk_free(&owner->mpi_owner_chunk);",
                      source, fixed = TRUE)
  status.pos <- max(gregexpr(
    "execution->status = NP_REGRESSION_GENERAL_LP_FIT_ERR_ROUTE;",
    source,
    fixed = TRUE
  )[[1L]])
  expect_gt(free.pos, 0L)
  expect_gt(status.pos, free.pos)
  expect_true(grepl("return R_NilValue;", source, fixed = TRUE))
})
