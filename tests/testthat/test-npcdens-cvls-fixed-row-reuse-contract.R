library(npRmpi)

locate_cvls_row_reuse_source <- function() {
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
  if (!length(hits)) NULL else hits[[1L]]
}

cvls_source_body <- function(lines, start_pattern, stop_pattern) {
  start <- grep(start_pattern, lines)
  stop <- grep(stop_pattern, lines)
  expect_gte(length(start), 1L)
  start <- tail(start, 1L)
  expect_gte(length(stop), 1L)
  stop <- stop[stop > start][1L]
  expect_true(is.finite(stop))
  paste(lines[start:(stop - 1L)], collapse = "\n")
}

.cvls_row_reuse_test_env <- new.env(parent = emptyenv())

.ensure_cvls_row_reuse_pool <- function() {
  if (!isTRUE(.cvls_row_reuse_test_env$started)) {
    npRmpi.init(nslaves = 1L, quiet = TRUE)
    .cvls_row_reuse_test_env$started <- TRUE
    withr::defer({
      if (isTRUE(.cvls_row_reuse_test_env$started)) {
        try(npRmpi.quit(force = TRUE), silent = TRUE)
        .cvls_row_reuse_test_env$started <- FALSE
      }
    }, envir = testthat::teardown_env())
  }
}

test_that("fixed and generalized-NN CVLS share one canonical LOO block engine", {
  src_file <- locate_cvls_row_reuse_source()
  skip_if(is.null(src_file), "source file src/jksum.c unavailable in this test context")
  lines <- readLines(src_file, warn = FALSE)
  source <- paste(lines, collapse = "\n")

  body <- cvls_source_body(
    lines,
    "^static int np_conditional_x_weight_block_stream_core_impl\\(",
    "^static int np_conditional_x_weight_block_stream_core\\("
  )

  expect_match(body, "BANDWIDTH_den_extern != BW_FIXED", fixed = TRUE)
  expect_match(body, "BANDWIDTH_den_extern != BW_GEN_NN", fixed = TRUE)
  expect_match(body, "const int use_bwctx =", fixed = TRUE)
  expect_match(
    body,
    "(BANDWIDTH_den_extern == BW_GEN_NN) ? matrix_bandwidth_x[l][i]",
    fixed = TRUE
  )
  expect_equal(lengths(regmatches(
    body,
    gregexpr("np_conditional_kernel_row_raw\\(", body, perl = TRUE)
  )), 1L)
  expect_false(grepl("np_glp_qr_drop_workspace_apply", body, fixed = TRUE))
  expect_false(grepl("self_weight = kw[eval_pos]", body, fixed = TRUE))
  markers <- c(
    "np_conditional_kernel_row_raw",
    "np_lp_full_row_workspace_solve",
    "np_lp_delete_denominator(rows_out[i][eval_idx], &den)"
  )
  positions <- vapply(markers, function(marker) {
    regexpr(marker, body, fixed = TRUE)[[1L]]
  }, integer(1L))
  expect_true(all(positions > 0L))
  expect_true(all(diff(positions) > 0L))
  expect_match(body, "(j == eval_pos) ? 0.0 : rows_out[i][orig_j]/den",
               fixed = TRUE)
  expect_false(grepl("np_conditional_x_weight_block_pair", source, fixed = TRUE))
})

test_that("conditional CVLS dispatch reaches only the canonical block owner", {
  src_file <- locate_cvls_row_reuse_source()
  skip_if(is.null(src_file), "source file src/jksum.c unavailable in this test context")
  lines <- readLines(src_file, warn = FALSE)

  density_body <- cvls_source_body(
    lines,
    "^int np_conditional_density_cvls_lp_stream\\(",
    "^static int np_conditional_distribution_cvls_lp_row_stream\\("
  )
  distribution_body <- cvls_source_body(
    lines,
    "^int np_conditional_distribution_cvls_lp_stream\\(",
    "^int np_kernel_estimate_density_categorical_leave_one_out_cv\\("
  )

  expect_equal(lengths(regmatches(
    density_body,
    gregexpr("np_conditional_x_weight_block_stream_core_suppress\\(", density_body, perl = TRUE)
  )), 2L)
  expect_equal(lengths(regmatches(
    density_body,
    gregexpr("np_conditional_x_weight_block_stream_core_impl\\(", density_body, perl = TRUE)
  )), 1L)
  expect_match(
    density_body,
    "((np_lp_engine_extern == NP_LP_ENGINE_SCALAR) ||",
    fixed = TRUE
  )
  expect_false(grepl(
    "np_conditional_x_weight_block_pair",
    distribution_body,
    fixed = TRUE
  ))

  shared_x_body <- cvls_source_body(
    lines,
    "^static int np_conditional_x_weight_block_stream_core_impl\\(",
    "^static int np_conditional_x_weight_block_stream_core\\("
  )
  expect_match(shared_x_body, "BANDWIDTH_den_extern != BW_FIXED", fixed = TRUE)
  expect_match(shared_x_body, "BANDWIDTH_den_extern != BW_GEN_NN", fixed = TRUE)
  expect_match(
    shared_x_body,
    "(BANDWIDTH_den_extern == BW_GEN_NN) ? matrix_bandwidth_x[l][i]",
    fixed = TRUE
  )

  shared_y_body <- cvls_source_body(
    lines,
    "^static int np_conditional_y_eval_block_stream_op_core\\(",
    "^#define NP_BOUNDED_CVLS_I1_GRID_POINTS"
  )
  expect_match(
    shared_y_body,
    "(BANDWIDTH_den_extern == BW_FIXED) ? 1 : block_rows",
    fixed = TRUE
  )

  expect_match(density_body, "block_id % iNum_Processors", fixed = TRUE)
  expect_match(density_body, "MPI_Allreduce(&local_fail, &any_fail", fixed = TRUE)
  expect_match(density_body, "MPI_Allreduce(MPI_IN_PLACE, block_terms", fixed = TRUE)
  expect_equal(lengths(regmatches(
    density_body,
    gregexpr("MPI_Allreduce\\(", density_body, perl = TRUE)
  )), 2L)
})

test_that("fixed CVLS row reuse preserves repeatable multivariate degree objectives", {
  skip_on_cran()
  .ensure_cvls_row_reuse_pool()
  set.seed(20260722)
  n <- 48L
  x <- data.frame(x1 = runif(n), x2 = runif(n))
  y <- data.frame(y1 = x$x1^2 - 0.3 * x$x2 + rnorm(n, sd = 0.08))
  bw <- npcdensbw(
    xdat = x,
    ydat = y,
    bws = c(0.38, 0.52, 0.44),
    bandwidth.compute = FALSE,
    bwmethod = "cv.ls",
    bwtype = "fixed",
    regtype = "lp",
    basis = "glp",
    degree = c(1L, 2L)
  )

  old_options <- options(np.tree = FALSE, np.categorical.compress = FALSE)
  on.exit(options(old_options), add = TRUE)
  first <- npRmpi:::.npcdensbw_eval_only(x, y, bw)
  second <- npRmpi:::.npcdensbw_eval_only(x, y, bw)

  expect_true(is.finite(first$objective))
  expect_identical(writeBin(first$objective, raw(), size = 8L),
                   writeBin(second$objective, raw(), size = 8L))
  expect_identical(first$num.feval, 1L)
  expect_identical(second$num.feval, 1L)
  expect_identical(as.integer(bw$degree.engine), c(1L, 2L))
})

test_that("generalized-NN CVLS row reuse preserves repeatable multivariate degree objectives", {
  skip_on_cran()
  .ensure_cvls_row_reuse_pool()
  set.seed(20260724)
  n <- 96L
  x <- data.frame(x1 = runif(n), x2 = runif(n))
  y <- data.frame(
    y1 = sin(2 * pi * x$x1) - 0.4 * x$x2^2 + rnorm(n, sd = 0.12)
  )
  bw <- npcdensbw(
    xdat = x,
    ydat = y,
    bws = c(31L, 29L, 37L),
    bandwidth.compute = FALSE,
    bwmethod = "cv.ls",
    bwtype = "generalized_nn",
    regtype = "lp",
    basis = "glp",
    degree = c(1L, 2L)
  )

  old_options <- options(np.tree = FALSE, np.categorical.compress = FALSE)
  on.exit(options(old_options), add = TRUE)
  first <- npRmpi:::.npcdensbw_eval_only(x, y, bw)
  second <- npRmpi:::.npcdensbw_eval_only(x, y, bw)

  expect_true(is.finite(first$objective))
  expect_identical(writeBin(first$objective, raw(), size = 8L),
                   writeBin(second$objective, raw(), size = 8L))
  expect_identical(first$num.feval, 1L)
  expect_identical(second$num.feval, 1L)
  expect_identical(as.integer(bw$degree.engine), c(1L, 2L))
})
