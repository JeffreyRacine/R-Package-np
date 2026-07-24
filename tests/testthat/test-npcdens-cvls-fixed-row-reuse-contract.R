library(npRmpi)

locate_fixed_cvls_row_reuse_source <- function() {
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

fixed_cvls_source_body <- function(lines, start_pattern, stop_pattern) {
  start <- grep(start_pattern, lines)
  stop <- grep(stop_pattern, lines)
  expect_gte(length(start), 1L)
  start <- tail(start, 1L)
  expect_gte(length(stop), 1L)
  stop <- stop[stop > start][1L]
  expect_true(is.finite(stop))
  paste(lines[start:(stop - 1L)], collapse = "\n")
}

test_that("fixed CVLS row reuse stays fixed-LP-only and preserves consumer order", {
  src_file <- locate_fixed_cvls_row_reuse_source()
  skip_if(is.null(src_file), "source file src/jksum.c unavailable in this test context")
  lines <- readLines(src_file, warn = FALSE)

  body <- fixed_cvls_source_body(
    lines,
    "^static int np_conditional_x_weight_block_pair_stream_core\\(",
    "^static int np_conditional_y_block_stream_op_core\\("
  )

  expect_match(body, "BANDWIDTH_den_extern != BW_FIXED", fixed = TRUE)
  expect_false(grepl("BW_GEN_NN", body, fixed = TRUE))
  expect_equal(lengths(regmatches(
    body,
    gregexpr("np_shadow_conditional_kernel_row_raw\\(", body, perl = TRUE)
  )), 1L)

  markers <- c(
    "self_weight = kw[eval_pos]",
    "kw[eval_pos] = 0.0",
    "np_glp_qr_drop_workspace_apply",
    "kw[eval_pos] = self_weight",
    "np_lp_full_row_workspace_solve"
  )
  positions <- vapply(markers, function(marker) {
    regexpr(marker, body, fixed = TRUE)[[1L]]
  }, integer(1L))
  expect_true(all(positions > 0L))
  expect_true(all(diff(positions) > 0L))
  expect_match(body, "suppress_nn_parallel", fixed = TRUE)
})

test_that("fixed CVLS row reuse leaves NN helpers and MPI reductions alone", {
  src_file <- locate_fixed_cvls_row_reuse_source()
  skip_if(is.null(src_file), "source file src/jksum.c unavailable in this test context")
  lines <- readLines(src_file, warn = FALSE)

  density_body <- fixed_cvls_source_body(
    lines,
    "^int np_conditional_density_cvls_lp_stream\\(",
    "^static int np_conditional_distribution_cvls_lp_row_stream\\("
  )

  expect_equal(lengths(regmatches(
    density_body,
    gregexpr("np_conditional_x_weight_block_pair_stream_core\\(", density_body, perl = TRUE)
  )), 1L)
  expect_match(density_body, "BANDWIDTH_den_extern == BW_FIXED", fixed = TRUE)
  expect_match(density_body, "int_ll_extern == LL_LP", fixed = TRUE)

  generalized_start <- regexpr(
    "if((BANDWIDTH_den_extern == BW_GEN_NN)",
    density_body,
    fixed = TRUE
  )[[1L]]
  pair_start <- regexpr(
    "np_conditional_x_weight_block_pair_stream_core(",
    density_body,
    fixed = TRUE
  )[[1L]]
  expect_true(generalized_start > 0L)
  expect_true(pair_start > generalized_start)
  generalized_block <- substr(density_body, generalized_start, pair_start - 1L)
  expect_false(grepl(
    "np_conditional_x_weight_block_pair_stream_core(",
    generalized_block,
    fixed = TRUE
  ))
  expect_equal(lengths(regmatches(
    generalized_block,
    gregexpr("np_conditional_x_weight_block_stream_core_impl\\(", generalized_block, perl = TRUE)
  )), 2L)

  shared_x_body <- fixed_cvls_source_body(
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
  src_file <- locate_fixed_cvls_row_reuse_source()
  skip_if(is.null(src_file), "source file src/jksum.c unavailable in this test context")
  lines <- readLines(src_file, warn = FALSE)

  distribution_body <- fixed_cvls_source_body(
    lines,
    "^int np_conditional_distribution_cvls_lp_stream\\(",
    "^static int np_shadow_conditional_build_y_matrix\\("
  )
  expect_false(grepl(
    "np_conditional_x_weight_block_pair_stream_core(",
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
