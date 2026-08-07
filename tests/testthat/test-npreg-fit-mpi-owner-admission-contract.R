npreg_fit_owner_source <- function() {
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

npreg_fit_owner_body <- function(source) {
  start <- regexpr(
    "static SEXP np_regression_general_lp_fit_execute(void *data)",
    source,
    fixed = TRUE
  )[[1L]]
  expect_gt(start, 0L)
  remainder <- substr(source, start, nchar(source))
  stop <- regexpr(
    "static int np_regression_general_lp_fit(",
    remainder,
    fixed = TRUE
  )[[1L]]
  expect_gt(stop, 0L)
  substr(source, start, start + stop - 2L)
}

test_that("general LP fit owner admits fixed, generalized NN, and adaptive NN", {
  source <- npreg_fit_owner_source()
  skip_if(is.null(source), "package C source unavailable")
  body <- npreg_fit_owner_body(source)

  expect_match(
    body,
    paste0(
      "((BANDWIDTH_reg == BW_FIXED) ||\n",
      "\t       (BANDWIDTH_reg == BW_GEN_NN) ||\n",
      "\t       (BANDWIDTH_reg == BW_ADAP_NN))"
    ),
    fixed = TRUE
  )
})

test_that("nearest-neighbor owner selects the canonical bandwidth shape", {
  source <- npreg_fit_owner_source()
  skip_if(is.null(source), "package C source unavailable")
  body <- npreg_fit_owner_body(source)

  expect_match(
    body,
    "(BANDWIDTH_reg == BW_GEN_NN) ?\n\t                   owner->matrix_bandwidth_eval :\n\t                   call->matrix_bandwidth,",
    fixed = TRUE
  )
  expect_match(
    body,
    "if(BANDWIDTH_reg != BW_ADAP_NN)",
    fixed = TRUE
  )
  expect_match(
    body,
    "(BANDWIDTH_reg == BW_GEN_NN) ?\n\t                  call->matrix_bandwidth[l][jj] :",
    fixed = TRUE
  )
  expect_match(
    body,
    paste0(
      "(BANDWIDTH_reg == BW_GEN_NN) ?\n",
      "\t                                     owner->matrix_bandwidth_eval :\n",
      "\t                                     call->matrix_bandwidth,"
    ),
    fixed = TRUE
  )
  expect_match(
    body,
    "call->do_merr && call->kernel_route == NULL",
    fixed = TRUE
  )
})

test_that("MPI LP owner clamps residual variance before signed normalization", {
  source <- npreg_fit_owner_source()
  skip_if(is.null(source), "package C source unavailable")
  body <- npreg_fit_owner_body(source)

  clamp <- regexpr("if(sigma2_owner <= 0.0) {", body, fixed = TRUE)[[1L]]
  division <- regexpr(
    "const double v_owner = sigma2_owner *",
    body,
    fixed = TRUE
  )[[1L]]
  expect_gt(clamp, 0L)
  expect_gt(division, clamp)
  expect_match(body, "out_owner[1] = 0.0;", fixed = TRUE)
})
