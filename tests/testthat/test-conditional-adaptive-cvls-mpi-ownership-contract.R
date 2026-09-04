conditional_adaptive_cvls_source <- function() {
  candidates <- c(
    test_path("..", "..", "src", "jksum.c"),
    test_path("..", "..", "..", "src", "jksum.c"),
    file.path(Sys.getenv("R_PACKAGE_SOURCE", ""), "src", "jksum.c"),
    file.path(getwd(), "src", "jksum.c")
  )
  hits <- unique(candidates[nzchar(candidates) & file.exists(candidates)])
  if (!length(hits)) return(NULL)
  paste(readLines(hits[[1L]], warn = FALSE), collapse = "\n")
}

conditional_adaptive_cvls_body <- function(source, start, stop) {
  first <- regexpr(start, source, fixed = TRUE)[[1L]]
  expect_gt(first, 0L)
  remainder <- substr(source, first, nchar(source))
  offset <- regexpr(stop, remainder, fixed = TRUE)[[1L]]
  expect_gt(offset, 0L)
  substr(remainder, 1L, offset - 1L)
}

test_that("adaptive conditional density has an isolated MPI owner", {
  source <- conditional_adaptive_cvls_source()
  skip_if(is.null(source), "package C source unavailable")
  body <- conditional_adaptive_cvls_body(
    source,
    "static int np_conditional_density_cvls_lp_adap_block_parallel_stream(",
    "#define NP_CDENS_ADAP_WIDTH3_NOT_BENEFICIAL"
  )

  expect_match(body, "const int ownership_stride = iNum_Processors;",
               fixed = TRUE)
  expect_match(body, "const int first_owned_block = my_rank;", fixed = TRUE)
  expect_match(body, "np_objective_outer_buffer_prepare(", fixed = TRUE)
  expect_match(body, "block_id[g] = first_block_id + g*ownership_stride;",
               fixed = TRUE)
  expect_match(body, "block_terms[block_id[g]] = quad[g] - 2.0*lin[g];",
               fixed = TRUE)
  expect_match(body, "np_objective_outer_buffer_finish(", fixed = TRUE)
  expect_match(body, '"NP_RMPI_INJECT_CDEN_CVLS_FAIL_RANK"', fixed = TRUE)
  expect_match(body, "np_blas_dgemm_tn_int(", fixed = TRUE)
  expect_false(grepl("num_obs*num_obs", body, fixed = TRUE))
  expect_false(grepl("num_obs * num_obs", body, fixed = TRUE))
  expect_false(grepl("alloc_matd(num_obs, num_obs)", body, fixed = TRUE))
})

test_that("adaptive density allocation fallback remains MPI-owned", {
  source <- conditional_adaptive_cvls_source()
  skip_if(is.null(source), "package C source unavailable")
  selector <- conditional_adaptive_cvls_body(
    source,
    "if(BANDWIDTH_den_extern == BW_ADAP_NN) {",
    "if(((int_TREE_X == NP_TREE_TRUE) || (int_TREE_Y == NP_TREE_TRUE))"
  )
  row_owner <- conditional_adaptive_cvls_body(
    source,
    "static int np_conditional_density_cvls_lp_row_parallel_stream(",
    "static int np_conditional_density_cvls_lp_adap_block_parallel_stream("
  )

  expect_match(
    selector,
    "np_conditional_density_cvls_lp_adap_block_parallel_stream(",
    fixed = TRUE
  )
  expect_match(
    selector,
    "np_conditional_density_cvls_lp_row_parallel_stream(",
    fixed = TRUE
  )
  expect_match(row_owner, "np_objective_outer_owned_rows(", fixed = TRUE)
  expect_match(row_owner, "contributions[i] = quad - 2.0*lin;",
               fixed = TRUE)
  expect_match(row_owner, "np_objective_outer_buffer_finish(", fixed = TRUE)
  expect_false(grepl(
    "return np_conditional_density_cvls_lp_row_stream(",
    selector,
    fixed = TRUE
  ))
})

test_that("adaptive block ownership is exactly once", {
  owned_blocks <- function(nblocks, nranks, rank) {
    width <- min(4L, (nblocks %/% nranks) + as.integer(nblocks %% nranks != 0L))
    width <- max(1L, width)
    first <- rank
    out <- integer()
    while (first < nblocks) {
      candidate <- first + seq.int(0L, width - 1L) * nranks
      out <- c(out, candidate[candidate < nblocks])
      first <- first + width * nranks
    }
    out
  }

  for (nblocks in 1:37) {
    for (nranks in 2:8) {
      claims <- unlist(lapply(0:(nranks - 1L), function(rank)
        owned_blocks(nblocks, nranks, rank)), use.names = FALSE)
      expect_identical(sort(claims), 0:(nblocks - 1L))
      expect_true(all(tabulate(claims + 1L, nbins = nblocks) == 1L))
    }
  }
})

test_that("adaptive conditional distribution has an isolated MPI owner", {
  source <- conditional_adaptive_cvls_source()
  skip_if(is.null(source), "package C source unavailable")
  body <- conditional_adaptive_cvls_body(
    source,
    "static int np_conditional_distribution_cvls_lp_adap_block_parallel_stream(",
    "static int np_conditional_distribution_cvls_lp_adap_block_stream("
  )

  expect_match(body, "const int ownership_stride = iNum_Processors;",
               fixed = TRUE)
  expect_match(body, "const int first_owned_block = my_rank;", fixed = TRUE)
  expect_match(body, "np_objective_outer_buffer_prepare(", fixed = TRUE)
  expect_match(body, "block_id[g] = first_block_id + g*ownership_stride;",
               fixed = TRUE)
  expect_match(body, "block_terms[block_id[g]] = block_sum[g];",
               fixed = TRUE)
  expect_match(body, "np_objective_outer_buffer_finish(", fixed = TRUE)
  expect_match(body, '"NP_RMPI_INJECT_CDIST_CVLS_FAIL_RANK"', fixed = TRUE)
  expect_match(body, "np_blas_dgemm_tn_int(", fixed = TRUE)
  expect_false(grepl("num_train*num_eval", body, fixed = TRUE))
  expect_false(grepl("num_train * num_eval", body, fixed = TRUE))
  expect_false(grepl("alloc_matd(num_train, num_eval)", body, fixed = TRUE))
})

test_that("adaptive distribution allocation fallback remains MPI-owned", {
  source <- conditional_adaptive_cvls_source()
  skip_if(is.null(source), "package C source unavailable")
  selector <- conditional_adaptive_cvls_body(
    source,
    "int np_conditional_distribution_cvls_lp_stream(",
    "static inline int np_density_cvml_beta_row_contribution("
  )
  row_owner <- conditional_adaptive_cvls_body(
    source,
    "static int np_conditional_distribution_cvls_lp_adap_row_parallel_stream(",
    "static int np_conditional_distribution_cvls_lp_adap_block_parallel_stream("
  )

  expect_match(
    selector,
    "np_conditional_distribution_cvls_lp_adap_block_parallel_stream(",
    fixed = TRUE
  )
  expect_match(
    selector,
    "np_conditional_distribution_cvls_lp_adap_row_parallel_stream(",
    fixed = TRUE
  )
  expect_match(row_owner, "np_objective_outer_owned_rows(", fixed = TRUE)
  expect_match(row_owner, "contributions[i] = row_loss;", fixed = TRUE)
  expect_match(row_owner, "np_objective_outer_buffer_finish(", fixed = TRUE)
})

test_that("exact adaptive conditional owners share rank-symmetric bounded finishing", {
  source <- conditional_adaptive_cvls_source()
  skip_if(is.null(source), "package C source unavailable")
  owners <- list(
    cvml = list(
      body = conditional_adaptive_cvls_body(
        source,
        "/* Exact ordinary adaptive-NN conditional-density CVML.",
        "/* Exact ordinary adaptive-NN conditional-density CVLS."
      ),
      hook = "NP_RMPI_INJECT_CDEN_ADAPTIVE_CVML_FAIL_RANK"
    ),
    cvls = list(
      body = conditional_adaptive_cvls_body(
        source,
        "/* Exact ordinary adaptive-NN conditional-density CVLS.",
        "int np_conditional_density_cvml_lp_stream("
      ),
      hook = "NP_RMPI_INJECT_CDEN_ADAPTIVE_CVLS_FAIL_RANK"
    ),
    cdist = list(
      body = conditional_adaptive_cvls_body(
        source,
        "/* Exact empirical-grid adaptive-NN conditional-distribution CV.",
        "/* Canonical generalized-NN empirical conditional-CDF row owner."
      ),
      hook = "NP_RMPI_INJECT_CDIST_ADAPTIVE_EXACT_FAIL_RANK"
    )
  )

  for (owner in owners) {
    body <- owner$body
    expect_match(
      body,
      "const int use_parallel_rows = np_objective_outer_rows_enabled(1);",
      fixed = TRUE
    )
    expect_match(body, "np_objective_outer_preflight_failed(", fixed = TRUE)
    expect_match(body, "np_objective_outer_buffer_prepare(", fixed = TRUE)
    expect_match(body, "np_objective_outer_owned_rows(", fixed = TRUE)
    expect_match(body, "np_objective_outer_buffer_finish(", fixed = TRUE)
    expect_match(body, owner$hook, fixed = TRUE)
    expect_false(grepl("num_obs*num_obs", body, fixed = TRUE))
    expect_false(grepl("num_train*num_train", body, fixed = TRUE))
    expect_false(grepl("alloc_matd(num_obs, num_obs)", body, fixed = TRUE))
    expect_false(grepl("alloc_matd(num_train, num_train)", body, fixed = TRUE))
  }
})

test_that("adaptive geometry preparation fails rank-symmetrically", {
  source <- conditional_adaptive_cvls_source()
  skip_if(is.null(source), "package C source unavailable")

  expect_match(
    source,
    "static int np_adaptive_geometry_preflight_failed(",
    fixed = TRUE
  )
  expect_match(
    source,
    "np_objective_outer_preflight_failed(\n    np_objective_outer_rows_enabled(1), local_fail);",
    fixed = TRUE
  )
  calls <- gregexpr(
    "np_adaptive_geometry_preflight_failed(", source, fixed = TRUE
  )[[1L]]
  # One helper definition plus five audited call sites: two regression,
  # one conditional-distribution, and two conditional-density owners.
  expect_identical(length(calls[calls > 0L]), 6L)
})
