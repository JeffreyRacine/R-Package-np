library(npRmpi)

locate_xrow_weighted_blas_source <- function() {
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

xrow_weighted_blas_source_body <- function(lines, start_pattern, stop_pattern) {
  start <- grep(start_pattern, lines)
  stop <- grep(stop_pattern, lines)
  expect_length(start, 1L)
  expect_gte(length(stop), 1L)
  stop <- stop[stop > start][1L]
  expect_true(is.finite(stop))
  paste(lines[start:(stop - 1L)], collapse = "\n")
}

test_that("conditional CVLS weighted BLAS gate is narrow and memory bounded", {
  src_file <- locate_xrow_weighted_blas_source()
  skip_if(is.null(src_file), "source file src/jksum.c unavailable in this test context")
  lines <- readLines(src_file, warn = FALSE)
  all_source <- paste(lines, collapse = "\n")
  gate <- xrow_weighted_blas_source_body(
    lines,
    "^static int np_apple_conditional_x_weighted_blas_profitable\\(",
    "^static void np_blas_dgemv_t_int\\("
  )

  expect_match(
    all_source,
    "#define NP_CONDITIONAL_X_WEIGHTED_BLAS_MIN_ROWS 2048",
    fixed = TRUE
  )
  expect_match(
    all_source,
    "#define NP_CONDITIONAL_X_WEIGHTED_BLAS_MAX_BYTES",
    fixed = TRUE
  )
  expect_match(
    all_source,
    "((size_t)64*(size_t)1024*(size_t)1024)",
    fixed = TRUE
  )
  expect_match(gate, "#if NP_ACCEL_GAUSS_COMPILED", fixed = TRUE)
  expect_match(gate, "!np_mseries_accelerate_enabled_cache", fixed = TRUE)
  expect_match(
    gate,
    "nrows < NP_CONDITIONAL_X_WEIGHTED_BLAS_MIN_ROWS",
    fixed = TRUE
  )
  expect_match(gate, "nterms < 2", fixed = TRUE)
  expect_match(gate, "basis_stride < nrows", fixed = TRUE)
  expect_match(gate, "SIZE_MAX/(size_t)nterms", fixed = TRUE)
  expect_match(
    gate,
    "NP_CONDITIONAL_X_WEIGHTED_BLAS_MAX_BYTES/sizeof(double)",
    fixed = TRUE
  )
  expect_match(
    gate,
    "np_apple_conditional_x_block_weighted_blas_profitable(",
    fixed = TRUE
  )
  expect_match(gate, "return (nterms >= 4) &&", fixed = TRUE)
  expect_match(gate, "(void)nrows;", fixed = TRUE)
  expect_match(gate, "(void)nterms;", fixed = TRUE)
  expect_match(gate, "(void)basis_stride;", fixed = TRUE)
  expect_false(grepl("num_train\\*num_train", gate))
  expect_false(grepl("nrows\\*nrows", gate))
})

test_that("canonical CVLS weighted BLAS preserves signed row algebra and fallback", {
  src_file <- locate_xrow_weighted_blas_source()
  skip_if(is.null(src_file), "source file src/jksum.c unavailable in this test context")
  lines <- readLines(src_file, warn = FALSE)
  body <- xrow_weighted_blas_source_body(
    lines,
    "^static int np_conditional_x_weight_block_stream_core_impl\\(",
    "^static int np_conditional_x_weight_block_stream_core_ctx\\("
  )
  compact <- gsub("[[:space:]]+", " ", body)

  expect_match(body, "BANDWIDTH_den_extern != BW_FIXED", fixed = TRUE)
  expect_match(body, "BANDWIDTH_den_extern != BW_GEN_NN", fixed = TRUE)
  expect_match(
    body,
    "np_apple_conditional_x_block_weighted_blas_profitable(",
    fixed = TRUE
  )
  expect_match(
    compact,
    "weighted_design = (double *)malloc(weighted_count*sizeof(double));",
    fixed = TRUE
  )
  expect_match(
    body,
    "use_weighted_blas = (weighted_design != NULL);",
    fixed = TRUE
  )
  expect_match(
    body,
    "weighted_row[j] = basis_row[j]*kw[j];",
    fixed = TRUE
  )
  expect_match(body, "F77_CALL(dgemm)", fixed = TRUE)
  expect_match(body, "F77_CALL(dgemv)", fixed = TRUE)
  expect_match(
    body,
    "rows_out[i][orig_j] = kw[j]*mean_row[j];",
    fixed = TRUE
  )
  expect_match(
    body,
    "np_lp_delete_denominator(rows_out[i][eval_idx], &den)",
    fixed = TRUE
  )
  expect_match(body, "if(weighted_design != NULL) free(weighted_design);", fixed = TRUE)
  expect_false(grepl("alloc_vecd\\(num_train\\)", body))
  expect_false(grepl("num_train\\*num_train", body))

  markers <- c(
    "np_conditional_kernel_row_raw",
    "weighted_row[j] = basis_row[j]*kw[j]",
    "F77_CALL(dgemm)",
    "np_lp_solve_workspace_solve_adjoint_ranked",
    "F77_CALL(dgemv)",
    "rows_out[i][orig_j] = kw[j]*mean_row[j]",
    "np_lp_delete_denominator(rows_out[i][eval_idx], &den)"
  )
  positions <- vapply(markers, function(marker) {
    regexpr(marker, body, fixed = TRUE)[[1L]]
  }, integer(1L))
  expect_true(all(positions > 0L))
  expect_true(all(diff(positions) > 0L))

  expect_match(
    compact,
    "if(use_weighted_blas){ const char trans_t",
    fixed = TRUE
  )
  expect_match(
    compact,
    "} else { for(l = 0; l < k; l++) for(j = 0; j < k; j++) solve_workspace.gram_source[l + j*k] = 0.0;",
    fixed = TRUE
  )
  expect_match(
    compact,
    "} else { for(j = 0; j < num_train; j++){ double zju = 0.0;",
    fixed = TRUE
  )
})

test_that("generalized-NN CVLS reaches the same weighted BLAS implementation", {
  src_file <- locate_xrow_weighted_blas_source()
  skip_if(is.null(src_file), "source file src/jksum.c unavailable in this test context")
  lines <- readLines(src_file, warn = FALSE)
  body <- xrow_weighted_blas_source_body(
    lines,
    "^static int np_conditional_x_weight_block_stream_core_impl\\(",
    "^static int np_conditional_x_weight_block_stream_core_ctx\\("
  )
  compact <- gsub("[[:space:]]+", " ", body)

  expect_match(body, "BANDWIDTH_den_extern != BW_GEN_NN", fixed = TRUE)
  expect_match(body, "BANDWIDTH_den_extern != BW_FIXED", fixed = TRUE)
  expect_match(
    body,
    "np_apple_conditional_x_block_weighted_blas_profitable(",
    fixed = TRUE
  )
  expect_match(
    compact,
    "weighted_design = (double *)malloc(weighted_count*sizeof(double));",
    fixed = TRUE
  )
  expect_match(
    body,
    "use_weighted_blas = (weighted_design != NULL);",
    fixed = TRUE
  )
  expect_match(
    body,
    "weighted_row[j] = basis_row[j]*kw[j];",
    fixed = TRUE
  )
  expect_match(body, "F77_CALL(dgemm)", fixed = TRUE)
  expect_match(body, "F77_CALL(dgemv)", fixed = TRUE)
  expect_match(
    body,
    "rows_out[i][orig_j] = kw[j]*mean_row[j];",
    fixed = TRUE
  )
  expect_match(
    body,
    "np_lp_delete_denominator(rows_out[i][eval_idx], &den)",
    fixed = TRUE
  )
  expect_match(body, "if(weighted_design != NULL) free(weighted_design);", fixed = TRUE)
  expect_false(grepl("alloc_vecd\\(num_train\\)", body))
  expect_false(grepl("num_train\\*num_train", body))

  markers <- c(
    "np_conditional_kernel_row_raw",
    "weighted_row[j] = basis_row[j]*kw[j]",
    "F77_CALL(dgemm)",
    "np_lp_solve_workspace_solve_adjoint_ranked",
    "F77_CALL(dgemv)",
    "rows_out[i][orig_j] = kw[j]*mean_row[j]",
    "np_lp_delete_denominator(rows_out[i][eval_idx], &den)"
  )
  positions <- vapply(markers, function(marker) {
    regexpr(marker, body, fixed = TRUE)[[1L]]
  }, integer(1L))
  expect_true(all(positions > 0L))
  expect_true(all(diff(positions) > 0L))

  expect_match(
    compact,
    "if(use_weighted_blas){ const char trans_t",
    fixed = TRUE
  )
  expect_match(
    compact,
    "} else { for(l = 0; l < k; l++) for(j = 0; j < k; j++) solve_workspace.gram_source[l + j*k] = 0.0;",
    fixed = TRUE
  )
  expect_match(
    compact,
    "} else { for(j = 0; j < num_train; j++){ double zju = 0.0;",
    fixed = TRUE
  )
})

test_that("shared conditional full-row blocks reuse bounded weighted BLAS", {
  src_file <- locate_xrow_weighted_blas_source()
  skip_if(is.null(src_file), "source file src/jksum.c unavailable in this test context")
  lines <- readLines(src_file, warn = FALSE)
  body <- xrow_weighted_blas_source_body(
    lines,
    "^static int np_conditional_x_weight_block_stream_core_impl\\(",
    "^static int np_conditional_x_weight_block_stream_core_ctx\\("
  )
  compact <- gsub("[[:space:]]+", " ", body)

  expect_match(
    body,
    "np_apple_conditional_x_block_weighted_blas_profitable(",
    fixed = TRUE
  )
  expect_match(
    compact,
    "weighted_design = (double *)malloc(weighted_count*sizeof(double));",
    fixed = TRUE
  )
  expect_match(body, "use_weighted_blas = (weighted_design != NULL);", fixed = TRUE)
  expect_match(body, "weighted_row[j] = basis_row[j]*kw[j];", fixed = TRUE)
  expect_match(body, "F77_CALL(dgemm)", fixed = TRUE)
  expect_match(body, "F77_CALL(dgemv)", fixed = TRUE)
  expect_match(body, "rows_out[i][orig_j] = kw[j]*mean_row[j];", fixed = TRUE)
  expect_match(body, "if(weighted_design != NULL) free(weighted_design);", fixed = TRUE)
  expect_false(grepl("num_train\\*num_train", body))

  expect_match(
    compact,
    "if(use_weighted_blas){ const char trans_t",
    fixed = TRUE
  )
  expect_match(
    compact,
    "} else { for(l = 0; l < k; l++) for(j = 0; j < k; j++) solve_workspace.gram_source[l + j*k] = 0.0;",
    fixed = TRUE
  )
  expect_match(
    compact,
    "} else { for(j = 0; j < num_train; j++){ double zju = 0.0;",
    fixed = TRUE
  )
  expect_match(
    body,
    "np_lp_delete_denominator(rows_out[i][eval_idx], &den)",
    fixed = TRUE
  )
})

test_that("conditional weighted BLAS preserves fixed-candidate objectives at boundary widths", {
  old_options <- options(
    np.messages = FALSE,
    np.tree = FALSE,
    np.largeh = FALSE
  )
  old_seed_exists <- exists(".Random.seed", envir = .GlobalEnv,
                            inherits = FALSE)
  if (old_seed_exists)
    old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit({
    options(old_options)
    if (old_seed_exists) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)

  set.seed(20260822L)
  seed_before <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  n <- 2048L
  i <- seq_len(n)
  xdat <- data.frame(x = (i - 0.5) / n)
  ydat <- data.frame(
    y = 0.35 * sin(2 * pi * xdat$x) +
      0.13 * cos(17 * pi * xdat$x) +
      0.08 * sin(i * sqrt(2))
  )
  cases <- data.frame(
    method = c("cv.ml", "cv.ml", "cv.ml", "cv.ls", "cv.ls"),
    degree = c(1L, 2L, 3L, 1L, 3L),
    stringsAsFactors = FALSE
  )

  evaluate <- function(state, accelerate) {
    options(np.macMseries.accelerate = accelerate)
    first <- second <- NULL
    expect_silent(first <- npRmpi:::.npcdensbw_eval_only(
      xdat, ydat, state, force.local = TRUE
    ))
    expect_silent(second <- npRmpi:::.npcdensbw_eval_only(
      xdat, ydat, state, force.local = TRUE
    ))
    expect_identical(first, second)
    expect_identical(names(first),
                     c("objective", "num.feval", "num.feval.fast",
                       "num.feval.guarded"))
    expect_true(is.finite(first$objective))
    expect_identical(first$num.feval, 1L)
    expect_identical(first$num.feval.fast, 0)
    expect_length(first$num.feval.guarded, 1L)
    expect_true(is.na(first$num.feval.guarded) ||
                (is.finite(first$num.feval.guarded) &&
                 first$num.feval.guarded >= 0))
    first
  }

  make_state <- function(method, degree) {
    reg_args <- list(
      bwmethod = method,
      bwscaling = FALSE,
      bwtype = "fixed",
      cxkertype = "gaussian",
      cxkerorder = 2L,
      cxkerbound = "none",
      uxkertype = "aitchisonaitken",
      oxkertype = "liracine",
      cykertype = "gaussian",
      cykerorder = 2L,
      cykerbound = "none",
      uykertype = "aitchisonaitken",
      oykertype = "liracine",
      regtype = "lp",
      pregtype = "Local-Polynomial",
      basis = "glp",
      degree = degree,
      bernstein.basis = FALSE,
      regtype.engine = "lp",
      basis.engine = "glp",
      degree.engine = degree,
      bernstein.basis.engine = FALSE,
      scale.factor.search.lower = 0.1,
      cvls.quadrature.grid = "hybrid",
      cvls.quadrature.extend.factor = 1,
      cvls.quadrature.points = c(100L, 50L),
      cvls.quadrature.ratios = c(0.20, 0.55, 0.25)
    )
    npRmpi:::.npcdensbw_build_conbandwidth(
      xdat = xdat,
      ydat = ydat,
      bws = c(0.22, 0.18),
      bandwidth.compute = FALSE,
      reg.args = reg_args
    )
  }

  for (case in seq_len(nrow(cases))) {
    method <- cases$method[[case]]
    degree <- cases$degree[[case]]
    state <- NULL
    expect_silent(state <- make_state(method, degree))

    expect_identical(state$method, method)
    expect_identical(state$type, "fixed")
    expect_identical(state$regtype.engine, "lp")
    expect_identical(state$degree.engine, degree)
    expect_identical(state$basis.engine, "glp")
    expect_false(state$bernstein.basis.engine)

    ordinary <- evaluate(state, FALSE)
    accelerated <- evaluate(state, TRUE)
    expect_identical(ordinary[names(ordinary) != "objective"],
                     accelerated[names(accelerated) != "objective"])
    scale <- max(1, abs(ordinary$objective), abs(accelerated$objective))
    error_bound <- 1024 * .Machine$double.eps * scale
    error <- abs(accelerated$objective - ordinary$objective)
    expect_true(
      error <= error_bound,
      info = sprintf(
        "%s width %d: abs.error=%.17g, bound=%.17g",
        method, degree + 1L, error, error_bound
      )
    )
  }

  expect_identical(
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE),
    seed_before
  )
})
