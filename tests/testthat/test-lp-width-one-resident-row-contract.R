locate_lp_row_source <- function() {
  candidates <- c(
    test_path("..", "..", "src", "jksum_lp_row.c"),
    test_path("..", "..", "..", "src", "jksum_lp_row.c"),
    file.path(Sys.getenv("R_PACKAGE_DIR", ""), "src", "jksum_lp_row.c"),
    file.path(Sys.getenv("R_PACKAGE_SOURCE", ""), "src", "jksum_lp_row.c"),
    file.path(getwd(), "src", "jksum_lp_row.c"),
    file.path(getwd(), "..", "src", "jksum_lp_row.c")
  )
  candidates <- unique(candidates[nzchar(candidates)])
  hits <- candidates[file.exists(candidates)]
  if (!length(hits))
    return(NULL)
  hits[[1L]]
}

test_that("fixed-objective width-one rows use an implicit scalar basis", {
  row_file <- locate_lp_row_source()
  skip_if(is.null(row_file), "source file src/jksum_lp_row.c unavailable")
  source <- paste(readLines(row_file, warn = FALSE), collapse = "\n")

  start <- regexpr(
    "static void np_lp_accumulate_dense_row_1(",
    source,
    fixed = TRUE
  )
  stop <- regexpr(
    "#define NP_LP_DEFINE_RESIDENT_WIDTH(WIDTH)",
    source,
    fixed = TRUE
  )
  expect_gt(start, 0L)
  expect_gt(stop, start)
  scalar <- substr(source, start, stop - 1L)

  expect_true(grepl("if(!ctx->use_tree)", scalar, fixed = TRUE))
  expect_true(grepl("fixed_rhs += weight*ctx->response[orig_ii]", scalar,
                    fixed = TRUE))
  expect_true(grepl("ctx->rhs[orig_ii] += weight*eval_response", scalar,
                    fixed = TRUE))
  expect_true(grepl("fixed_moment += weight", scalar, fixed = TRUE))
  expect_true(grepl("ctx->moments[orig_ii] += weight", scalar,
                    fixed = TRUE))
  expect_false(grepl("support", scalar, fixed = TRUE))
  expect_false(grepl("ctx->basis", scalar, fixed = TRUE))
  expect_false(grepl("F77_", scalar, fixed = TRUE))
  expect_false(grepl("LL_LC", scalar, fixed = TRUE))
  expect_false(grepl("NP_LP_DEFINE_RESIDENT_WIDTH(1)", source, fixed = TRUE))
  expect_true(grepl("NP_LP_DEFINE_RESIDENT_WIDTH(2)", source, fixed = TRUE))
})

test_that("sparse-tree width-one pairs use the same implicit scalar algebra", {
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
  skip_if(!length(hits), "source file src/jksum.c unavailable")
  source <- paste(readLines(hits[[1L]], warn = FALSE), collapse = "\n")

  start <- regexpr("static inline void np_lp_accumulate_pair(", source,
                   fixed = TRUE)
  stop <- regexpr(
    "static inline void np_lp_accumulate_sparse_pair_wide_resident(",
    source,
    fixed = TRUE
  )
  expect_gt(start, 0L)
  expect_gt(stop, start)
  pair <- substr(source, start, stop - 1L)

  scalar_start <- regexpr("if(nterms == 1)", pair, fixed = TRUE)
  wider_start <- regexpr("if(nterms == 2)", pair, fixed = TRUE)
  expect_gt(scalar_start, 0L)
  expect_gt(wider_start, scalar_start)
  scalar <- substr(pair, scalar_start, wider_start - 1L)
  expect_true(grepl("tj[0] += w*yi", scalar, fixed = TRUE))
  expect_true(grepl("ti[0] += w*eval_ybasis[0]", scalar, fixed = TRUE))
  expect_true(grepl("sj[0] += w", scalar, fixed = TRUE))
  expect_true(grepl("si[0] += w", scalar, fixed = TRUE))
  expect_false(grepl("support", scalar, fixed = TRUE))
  expect_false(grepl("(^|[^[:alnum:]_])basis\\[", scalar, perl = TRUE))
  expect_false(grepl("F77_", scalar, fixed = TRUE))
})

test_that("raw and Bernstein degree-zero objectives retain one scalar result", {
  package_owner <- environmentName(environment(npregbw))
  skip_if(
    identical(package_owner, "npRmpi"),
    "npRmpi numerical equality is covered by active-pool parity tests"
  )
  eval_only <- get(
    ".npregbw_eval_only",
    envir = asNamespace(package_owner),
    inherits = FALSE
  )

  set.seed(20260725)
  n <- 96L
  x <- data.frame(x1 = runif(n), x2 = runif(n))
  y <- sin(2 * pi * x$x1) * cos(2 * pi * x$x2) +
    rnorm(n, sd = 0.2)

  make_bw <- function(bernstein) {
    npregbw(
      xdat = x,
      ydat = y,
      bws = c(0.22, 0.24),
      regtype = "lp",
      basis = "glp",
      degree = c(0L, 0L),
      bernstein.basis = bernstein,
      bwmethod = "cv.ls",
      bandwidth.compute = FALSE
    )
  }

  raw <- eval_only(x, y, make_bw(FALSE))$objective
  bernstein <- eval_only(x, y, make_bw(TRUE))$objective
  expect_identical(raw, bernstein)
})

test_that("scalar fixed objectives enter the resident LP0 owner without a basis", {
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
  skip_if(!length(hits), "source file src/jksum.c unavailable")
  source <- paste(readLines(hits[[1L]], warn = FALSE), collapse = "\n")

  sparse_start <- regexpr(
    "np_lp_fixed_tree_sparse_accumulate(",
    source,
    fixed = TRUE
  )
  fixed_start <- regexpr(
    "static NPRegCvLpResult np_regression_cv_lp_basis_fixed(",
    source,
    fixed = TRUE
  )
  objective_start <- regexpr(
    "double np_kernel_estimate_regression_categorical_ls_aic(",
    source,
    fixed = TRUE
  )
  expect_gt(sparse_start, 0L)
  expect_gt(fixed_start, sparse_start)
  expect_gt(objective_start, fixed_start)
  sparse_owner <- substr(source, sparse_start, fixed_start - 1L)
  fixed_owner <- substr(source, fixed_start, objective_start - 1L)

  expect_true(grepl(
    "if((nterms <= 0) || ((nterms > 1) && (basis == NULL)))",
    fixed_owner,
    fixed = TRUE
  ))
  expect_true(grepl("eval_basis[0] = 1.0", fixed_owner, fixed = TRUE))
  expect_true(grepl("eval_ybasis[0] = yj", sparse_owner, fixed = TRUE))
  expect_true(grepl("eval_outer[0] = 1.0", sparse_owner, fixed = TRUE))
  expect_equal(
    lengths(regmatches(
      sparse_owner,
      gregexpr("np_lp_accumulate_pair_width_one(", sparse_owner, fixed = TRUE)
    )),
    2L
  )
  expect_true(grepl("scalar_moment = moments[j]", sparse_owner, fixed = TRUE))
  expect_true(grepl("scalar_rhs = rhs[j]", sparse_owner, fixed = TRUE))
  expect_equal(
    lengths(regmatches(
      sparse_owner,
      gregexpr("moments[j] = scalar_moment", sparse_owner, fixed = TRUE)
    )),
    2L
  )
  expect_equal(
    lengths(regmatches(
      sparse_owner,
      gregexpr("rhs[j] = scalar_rhs", sparse_owner, fixed = TRUE)
    )),
    2L
  )

  objective_end <- regexpr(
    "int kernel_estimate_regression_categorical_tree_np(",
    source,
    fixed = TRUE
  )
  expect_gt(objective_end, objective_start)
  objective <- substr(source, objective_start, objective_end - 1L)
  scalar_start <- regexpr(
    "if((lp_engine == NP_LP_ENGINE_SCALAR)",
    objective,
    fixed = TRUE
  )
  legacy_start <- regexpr(
    "if(lp_engine == NP_LP_ENGINE_SCALAR) {",
    objective,
    fixed = TRUE
  )
  expect_gt(scalar_start, 0L)
  expect_gt(legacy_start, scalar_start)
  scalar_owner <- substr(objective, scalar_start, legacy_start - 1L)

  expect_true(grepl(
    "np_reg_cv_scalar_use_resident_fixed(bwm,",
    scalar_owner,
    fixed = TRUE
  ))
  expect_true(grepl(
    "((!ks_tree_use) || (num_reg_continuous <= 2))",
    source,
    fixed = TRUE
  ))
  expect_true(grepl(
    "((bwm == RBWM_CVKS) && ks_tree_use &&",
    source,
    fixed = TRUE
  ))
  expect_true(grepl(
    "np_regression_cv_lp_basis_fixed(bwm",
    scalar_owner,
    fixed = TRUE
  ))
  expect_true(grepl(
    "                                      1,\n                                      NULL,",
    scalar_owner,
    fixed = TRUE
  ))
  expect_false(grepl("F77_", scalar_owner, fixed = TRUE))
})

test_that("conditional width-one blocks share one raw row without LAPACK", {
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
  skip_if(!length(hits), "source file src/jksum.c unavailable")
  source <- paste(readLines(hits[[1L]], warn = FALSE), collapse = "\n")

  lines <- strsplit(source, "\n", fixed = TRUE)[[1L]]
  start <- grep(
    "^static int np_conditional_x_weight_block_stream_core_impl\\(",
    lines
  )
  stop <- grep("^static int np_conditional_x_weight_block_stream_core\\(", lines)
  expect_length(start, 1L)
  expect_length(stop, 1L)
  expect_lt(start, stop)
  body <- paste(lines[start:(stop - 1L)], collapse = "\n")

  expect_true(grepl(
    "if(drop_eval_self && (lp_engine == NP_LP_ENGINE_SCALAR))",
    body,
    fixed = TRUE
  ))
  expect_true(grepl("kw[eval_pos] = 0.0;", body, fixed = TRUE))
  expect_true(grepl("if(lp_engine == NP_LP_ENGINE_SCALAR){", body,
                    fixed = TRUE))
  expect_true(grepl("rows_out[i][orig_j] = kw[j]/row_sum;", body,
                    fixed = TRUE))
  expect_false(grepl("np_conditional_x_weight_block_pair", source,
                     fixed = TRUE))
  expect_false(grepl("np_glp_qr_drop_workspace", body, fixed = TRUE))
})

test_that("conditional scalar streams select their canonical topology owner", {
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
  skip_if(!length(hits), "source file src/jksum.c unavailable")
  source <- paste(readLines(hits[[1L]], warn = FALSE), collapse = "\n")

  support_start <- regexpr(
    "int np_conditional_lp_stream_engine_supported(void){",
    source,
    fixed = TRUE
  )
  cvml_start <- regexpr(
    "int np_conditional_density_cvml_lp_stream(",
    source,
    fixed = TRUE
  )
  expect_gt(support_start, 0L)
  expect_gt(cvml_start, support_start)
  support <- substr(source, support_start, cvml_start - 1L)

  expect_true(grepl(
    "np_lp_engine_extern == NP_LP_ENGINE_SCALAR",
    support,
    fixed = TRUE
  ))
  expect_true(grepl("num_reg_continuous_extern > 0", support, fixed = TRUE))
  expect_true(grepl(
    "BANDWIDTH_den_extern != BW_ADAP_NN",
    support,
    fixed = TRUE
  ))
  cvml_support_start <- regexpr(
    "int np_conditional_density_cvml_stream_engine_supported(void){",
    support,
    fixed = TRUE
  )
  expect_gt(cvml_support_start, 0L)
  cvml_support <- substr(
    support,
    cvml_support_start,
    nchar(support)
  )
  expect_true(grepl(
    "num_reg_continuous_extern > 1",
    cvml_support,
    fixed = TRUE
  ))
  expect_true(grepl(
    "np_conditional_density_cvml_scalar_sparse_batch_supported()",
    cvml_support,
    fixed = TRUE
  ))
  expect_true(grepl(
    "never public regtype spelling or measured run time",
    cvml_support,
    fixed = TRUE
  ))
  sparse_support_start <- regexpr(
    "np_conditional_density_cvml_scalar_sparse_batch_supported(void){",
    support,
    fixed = TRUE
  )
  expect_gt(sparse_support_start, 0L)
  sparse_support <- substr(
    support,
    sparse_support_start,
    cvml_support_start - 1L
  )
  expect_true(grepl(
    "BANDWIDTH_den_extern == BW_FIXED",
    sparse_support,
    fixed = TRUE
  ))
  expect_true(grepl(
    "num_reg_continuous_extern >= 2",
    sparse_support,
    fixed = TRUE
  ))
  expect_true(grepl(
    "num_reg_continuous_extern <= 3",
    sparse_support,
    fixed = TRUE
  ))
  expect_true(grepl("int_TREE_X == NP_TREE_TRUE", sparse_support, fixed = TRUE))
  expect_true(grepl(
    "Width one has two canonical serial microkernels",
    cvml_support,
    fixed = TRUE
  ))
  expect_true(grepl(
    "if(!np_conditional_density_cvml_stream_engine_supported())",
    source,
    fixed = TRUE
  ))
  expect_true(grepl(
    "np_conditional_x_weight_block_stream_core(",
    source,
    fixed = TRUE
  ))
})
