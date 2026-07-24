library(npRmpi)

locate_yconv_supertile_source <- function() {
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

yconv_source_body <- function(lines, start_pattern, stop_pattern) {
  start <- grep(start_pattern, lines)
  stop <- grep(stop_pattern, lines)
  expect_gte(length(start), 1L)
  start <- tail(start, 1L)
  expect_gte(length(stop), 1L)
  stop <- stop[stop > start][1L]
  expect_true(is.finite(stop))
  paste(lines[start:(stop - 1L)], collapse = "\n")
}

test_that("MPI CVLS Y convolution supertile is memory bounded and isolated", {
  src_file <- locate_yconv_supertile_source()
  skip_if(is.null(src_file), "source file src/jksum.c unavailable in this test context")
  lines <- readLines(src_file, warn = FALSE)
  body <- yconv_source_body(
    lines,
    "^static int np_conditional_density_cvls_lp_supertile2_stream\\(",
    "^int np_conditional_density_cvls_lp_stream\\("
  )

  expect_match(body, "nblocks <= 1", fixed = TRUE)
  expect_match(body, "int_ll_extern != LL_LP", fixed = TRUE)
  expect_match(body, "num_reg_continuous_extern <= 0", fixed = TRUE)
  expect_match(body, "vector_glp_degree_extern == NULL", fixed = TRUE)
  expect_match(body, "BANDWIDTH_den_extern != BW_FIXED", fixed = TRUE)
  expect_match(body, "BANDWIDTH_den_extern != BW_GEN_NN", fixed = TRUE)
  expect_match(body, "int_TREE_X == NP_TREE_TRUE", fixed = TRUE)
  expect_match(body, "int_TREE_Y == NP_TREE_TRUE", fixed = TRUE)
  expect_false(grepl("bernstein", body, ignore.case = TRUE))

  expect_equal(lengths(regmatches(
    body,
    gregexpr("alloc_tmatd\\(num_obs, block_size\\)", body, perl = TRUE)
  )), 4L)
  expect_equal(lengths(regmatches(
    body,
    gregexpr("alloc_vecd\\(block_size\\*block_size\\)", body, perl = TRUE)
  )), 1L)
  expect_false(grepl("num_obs\\*num_obs", body))
  expect_false(grepl("num_obs \\* num_obs", body, fixed = TRUE))
})

test_that("MPI CVLS supertile retains rank ownership and block-order reduction", {
  src_file <- locate_yconv_supertile_source()
  skip_if(is.null(src_file), "source file src/jksum.c unavailable in this test context")
  lines <- readLines(src_file, warn = FALSE)
  body <- yconv_source_body(
    lines,
    "^static int np_conditional_density_cvls_lp_supertile2_stream\\(",
    "^int np_conditional_density_cvls_lp_stream\\("
  )

  expect_match(
    body,
    "use_parallel_blocks ? iNum_Processors : 1",
    fixed = TRUE
  )
  expect_match(body, "use_parallel_blocks ? my_rank : 0", fixed = TRUE)
  expect_match(
    body,
    "first_block_id += 2*ownership_stride",
    fixed = TRUE
  )
  expect_match(
    body,
    "second_block_id = first_block_id + ownership_stride",
    fixed = TRUE
  )
  expect_match(
    body,
    "block_terms[first_block_id]",
    fixed = TRUE
  )
  expect_match(
    body,
    "block_terms[second_block_id]",
    fixed = TRUE
  )

  markers <- c(
    "lin_first += np_blas_ddot_int",
    "lin_second += np_blas_ddot_int",
    "for(j0 = 0; j0 < num_obs; j0 += block_size)",
    "quad_first += aij*quad_cross",
    "quad_second += aij*quad_cross",
    "block_terms[first_block_id]",
    "block_terms[second_block_id]",
    "MPI_Allreduce(&local_fail, &any_fail",
    "MPI_Allreduce(MPI_IN_PLACE",
    "for(ii = 0; ii < nblocks; ii++)"
  )
  positions <- vapply(markers, function(marker) {
    regexpr(marker, body, fixed = TRUE)[[1L]]
  }, integer(1L))
  expect_true(all(positions > 0L))
  expect_true(all(diff(positions) > 0L))
  expect_equal(lengths(regmatches(
    body,
    gregexpr("MPI_Allreduce\\(", body, perl = TRUE)
  )), 2L)
})

test_that("MPI CVLS supertile dispatch leaves no-gain and excluded routes intact", {
  src_file <- locate_yconv_supertile_source()
  skip_if(is.null(src_file), "source file src/jksum.c unavailable in this test context")
  lines <- readLines(src_file, warn = FALSE)
  body <- yconv_source_body(
    lines,
    "^int np_conditional_density_cvls_lp_stream\\(",
    "^static int np_conditional_distribution_cvls_lp_row_stream\\("
  )

  dispatch <- paste(
    "if((int_ll_extern == LL_LP) &&",
    "(num_reg_continuous_extern > 0) &&",
    "(vector_glp_degree_extern != NULL) &&",
    "(int_TREE_X != NP_TREE_TRUE) &&",
    "(int_TREE_Y != NP_TREE_TRUE) &&",
    "(nblocks > 1) &&",
    "((!use_parallel_blocks) || (nblocks > iNum_Processors)))"
  )
  expect_match(gsub("[[:space:]]+", " ", body), dispatch, fixed = TRUE)
  expect_equal(lengths(regmatches(
    body,
    gregexpr(
      "np_conditional_density_cvls_lp_supertile2_stream\\(",
      body,
      perl = TRUE
    )
  )), 1L)

  adaptive_pos <- regexpr(
    "if(BANDWIDTH_den_extern == BW_ADAP_NN)",
    body,
    fixed = TRUE
  )[[1L]]
  tree_pos <- regexpr(
    "return np_conditional_density_cvls_lp_row_stream",
    body,
    fixed = TRUE
  )[[1L]]
  supertile_pos <- regexpr(
    "return np_conditional_density_cvls_lp_supertile2_stream",
    body,
    fixed = TRUE
  )[[1L]]
  allocation_pos <- regexpr(
    "xblock = alloc_tmatd",
    body,
    fixed = TRUE
  )[[1L]]
  expect_true(
    adaptive_pos > 0L &&
      tree_pos > adaptive_pos &&
      supertile_pos > tree_pos &&
      allocation_pos > supertile_pos
  )
})

test_that("target-active one-slave raw and Bernstein supertiles repeat", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(), add = TRUE)

  set.seed(20260724)
  n <- 1025L
  x <- data.frame(x1 = runif(n), x2 = runif(n))
  y <- data.frame(
    y1 = sin(2 * pi * x$x1) - 0.3 * x$x2^2 + rnorm(n, sd = 0.18)
  )
  old_options <- options(np.tree = FALSE, np.categorical.compress = FALSE)
  on.exit(options(old_options), add = TRUE)

  evaluate <- function(bwtype, bernstein) {
    bws <- if (identical(bwtype, "fixed")) {
      c(0.42, 0.38, 0.46)
    } else {
      c(307L, 281L, 337L)
    }
    bw <- npcdensbw(
      xdat = x,
      ydat = y,
      bws = bws,
      bandwidth.compute = FALSE,
      bwmethod = "cv.ls",
      bwtype = bwtype,
      regtype = "lp",
      basis = "glp",
      degree = c(2L, 1L),
      degree.select = "manual",
      bernstein.basis = bernstein
    )
    first <- npRmpi:::.npcdensbw_eval_only(x, y, bw)
    second <- npRmpi:::.npcdensbw_eval_only(x, y, bw)
    expect_identical(
      writeBin(first$objective, raw(), size = 8L),
      writeBin(second$objective, raw(), size = 8L)
    )
    expect_identical(first$num.feval, 1L)
    first$objective
  }

  for (bwtype in c("fixed", "generalized_nn")) {
    raw <- evaluate(bwtype, FALSE)
    bernstein <- evaluate(bwtype, TRUE)
    expect_equal(raw, bernstein, tolerance = 1e-11)
  }
})
