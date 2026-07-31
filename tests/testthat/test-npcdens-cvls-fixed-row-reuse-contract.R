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

test_that("fixed and generalized-NN CVLS derive LOO from their full rows", {
  src_file <- locate_cvls_row_reuse_source()
  skip_if(is.null(src_file), "source file src/jksum.c unavailable in this test context")
  lines <- readLines(src_file, warn = FALSE)

  fixed_body <- cvls_source_body(
    lines,
    "^static int (NP_NOINLINE )?np_conditional_x_weight_block_pair_stream_core\\(",
    "^static int (NP_NOINLINE )?(NP_HOT_ALIGN )?np_conditional_x_weight_block_pair_gnn_stream_core\\("
  )
  gnn_body <- cvls_source_body(
    lines,
    "^static int (NP_NOINLINE )?(NP_HOT_ALIGN )?np_conditional_x_weight_block_pair_gnn_stream_core\\(",
    "^/\\* NP_CONDITIONAL_X_WEIGHT_BLOCK_PAIR_GNN_STREAM_CORE_END \\*/$"
  )

  expect_match(fixed_body, "BANDWIDTH_den_extern != BW_FIXED", fixed = TRUE)
  expect_false(grepl("BW_GEN_NN", fixed_body, fixed = TRUE))
  expect_match(gnn_body, "BANDWIDTH_den_extern != BW_GEN_NN", fixed = TRUE)
  expect_false(grepl("BANDWIDTH_den_extern != BW_FIXED", gnn_body, fixed = TRUE))
  expect_match(
    gnn_body,
    "matrix_bandwidth_eval_one[l][0] = matrix_bandwidth_x[l][i]",
    fixed = TRUE
  )
  for (body in list(fixed_body, gnn_body)) {
    expect_equal(lengths(regmatches(
      body,
      gregexpr("np_shadow_conditional_kernel_row_raw\\(", body, perl = TRUE)
    )), 1L)
    expect_match(body, "suppress_nn_parallel", fixed = TRUE)
  }

  for (body in list(fixed_body, gnn_body)) {
    expect_false(grepl("np_glp_qr_drop_workspace_apply", body, fixed = TRUE))
    expect_false(grepl("self_weight = kw[eval_pos]", body, fixed = TRUE))
    expect_match(
      body,
      "np_lp_delete_denominator(full_rows_out[i][eval_idx], &den)",
      fixed = TRUE
    )
  }
  fixed_markers <- c(
    "np_lp_full_row_workspace_solve",
    "np_lp_delete_denominator(full_rows_out[i][eval_idx], &den)",
    "loo_rows_out[i][orig_j] ="
  )
  fixed_positions <- vapply(fixed_markers, function(marker) {
    regexpr(marker, fixed_body, fixed = TRUE)[[1L]]
  }, integer(1L))
  expect_true(all(fixed_positions > 0L))
  expect_true(all(diff(fixed_positions) > 0L))

  gnn_markers <- c(
    "matrix_bandwidth_eval_one[l][0] = matrix_bandwidth_x[l][i]",
    "np_shadow_conditional_kernel_row_raw",
    "np_lp_full_row_workspace_solve",
    "np_lp_delete_denominator(full_rows_out[i][eval_idx], &den)",
    "loo_rows_out[i][j] ="
  )
  gnn_positions <- vapply(gnn_markers, function(marker) {
    regexpr(marker, gnn_body, fixed = TRUE)[[1L]]
  }, integer(1L))
  expect_true(all(gnn_positions > 0L))
  expect_true(all(diff(gnn_positions) > 0L))
})

test_that("CVLS row reuse selects scalar, fixed, and generalized-NN siblings centrally", {
  src_file <- locate_cvls_row_reuse_source()
  skip_if(is.null(src_file), "source file src/jksum.c unavailable in this test context")
  lines <- readLines(src_file, warn = FALSE)

  density_body <- cvls_source_body(
    lines,
    "^int np_conditional_density_cvls_lp_stream\\(",
    "^static int np_conditional_distribution_cvls_lp_row_stream\\("
  )

  expect_equal(lengths(regmatches(
    density_body,
    gregexpr("np_conditional_x_weight_block_pair_selected_stream_core\\(", density_body, perl = TRUE)
  )), 2L)
  expect_equal(lengths(regmatches(
    density_body,
    gregexpr("np_conditional_x_weight_block_pair_stream_core\\(", density_body, perl = TRUE)
  )), 0L)
  expect_equal(lengths(regmatches(
    density_body,
    gregexpr("np_conditional_x_weight_block_pair_gnn_stream_core\\(", density_body, perl = TRUE)
  )), 0L)
  expect_match(
    density_body,
    "((np_lp_engine_extern == NP_LP_ENGINE_SCALAR) ||",
    fixed = TRUE
  )

  selected_body <- cvls_source_body(
    lines,
    "^static int np_conditional_x_weight_block_pair_selected_stream_core\\(",
    "^/\\* NP_CONDITIONAL_X_WEIGHT_BLOCK_PAIR_GNN_STREAM_CORE_END \\*/$"
  )
  markers <- c(
    "if(np_lp_engine_extern == NP_LP_ENGINE_SCALAR)",
    "return np_conditional_x_weight_block_pair_scalar_stream_core",
    "if(np_lp_engine_extern != NP_LP_ENGINE_GENERAL)",
    "if(BANDWIDTH_den_extern == BW_FIXED)",
    "return np_conditional_x_weight_block_pair_stream_core",
    "if(BANDWIDTH_den_extern == BW_GEN_NN)",
    "return np_conditional_x_weight_block_pair_gnn_stream_core"
  )
  positions <- vapply(markers, function(marker) {
    regexpr(marker, selected_body, fixed = TRUE)[[1L]]
  }, integer(1L))
  expect_true(all(positions > 0L))
  expect_true(all(diff(positions) > 0L))
  expect_false(grepl(
    "np_conditional_x_weight_block_stream_core_impl(",
    selected_body,
    fixed = TRUE
  ))

  shared_x_body <- cvls_source_body(
    lines,
    "^static int np_conditional_x_weight_block_stream_core_impl\\(",
    "^static int np_conditional_x_weight_block_stream_core\\("
  )
  expect_match(shared_x_body, "BANDWIDTH_den_extern != BW_GEN_NN", fixed = TRUE)
  expect_match(
    shared_x_body,
    "(BANDWIDTH_den_extern == BW_GEN_NN) ? matrix_bandwidth_x[l][i]",
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

test_that("fixed CVLS row reuse is isolated from npcdist", {
  src_file <- locate_cvls_row_reuse_source()
  skip_if(is.null(src_file), "source file src/jksum.c unavailable in this test context")
  lines <- readLines(src_file, warn = FALSE)

  distribution_body <- cvls_source_body(
    lines,
    "^int np_conditional_distribution_cvls_lp_stream\\(",
    "^int np_kernel_estimate_density_categorical_leave_one_out_cv\\("
  )
  expect_false(grepl(
    "np_conditional_x_weight_block_pair_stream_core(",
    distribution_body,
    fixed = TRUE
  ))
  expect_false(grepl(
    "np_conditional_x_weight_block_pair_gnn_stream_core(",
    distribution_body,
    fixed = TRUE
  ))
})

test_that("fixed CVLS row reuse preserves repeatable multivariate degree objectives", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(), add = TRUE)

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
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(), add = TRUE)

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
