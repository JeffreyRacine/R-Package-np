locate_conditional_xrow_source <- function() {
  roots <- c(
    test_path("..", ".."),
    test_path("..", "..", ".."),
    Sys.getenv("R_PACKAGE_DIR", ""),
    Sys.getenv("R_PACKAGE_SOURCE", ""),
    getwd(),
    file.path(getwd(), "..")
  )
  roots <- unique(roots[nzchar(roots)])
  roots <- roots[file.exists(file.path(roots, "src", "jksum.c"))]
  if (!length(roots)) return(NULL)
  file.path(roots[[1L]], "src", "jksum.c")
}

test_that("adaptive conditional X reciprocals remain a lazy optional sidecar", {
  path <- locate_conditional_xrow_source()
  skip_if(is.null(path), "package C sources unavailable in this test context")
  source <- paste(readLines(path, warn = FALSE), collapse = "\n")

  expect_match(
    source,
    "NPConditionalXRowReciprocalCache *reciprocal_cache;",
    fixed = TRUE
  )
  expect_match(
    source,
    "if(ctx->reciprocal_cache != NULL) free(ctx->reciprocal_cache);",
    fixed = TRUE
  )
  expect_match(
    source,
    "(np_glp_cv_cache.nterms < 4)",
    fixed = TRUE
  )
  expect_match(source, "(num_train <= 0) || (ndim < 2)", fixed = TRUE)
  expect_match(
    source,
    "(num_reg_unordered_extern != 0)",
    fixed = TRUE
  )
  expect_match(source, "int_cxker_bound_extern", fixed = TRUE)
  expect_match(
    source,
    "(!np_mseries_accelerate_enabled_cache)",
    fixed = TRUE
  )
  expect_match(
    source,
    "reciprocal_count = ((size_t)ndim + 1)*(size_t)num_train;",
    fixed = TRUE
  )

  prepare_start <- regexpr(
    "static int np_conditional_xrow_ctx_prepare",
    source,
    fixed = TRUE
  )
  row_start <- regexpr(
    "static int np_conditional_xrow_from_ctx_impl",
    source,
    fixed = TRUE
  )
  expect_gt(prepare_start, 0L)
  expect_gt(row_start, prepare_start)
  prepare <- substr(source, prepare_start, row_start - 1L)
  expect_false(grepl(
    "np_conditional_xrow_reciprocal_cache_try(ctx)",
    prepare,
    fixed = TRUE
  ))

  row_end <- regexpr(
    "static int np_conditional_xrow_from_ctx(",
    source,
    fixed = TRUE
  )
  expect_gt(row_end, row_start)
  row <- substr(source, row_start, row_end - 1L)
  expect_match(row, "(eval_idx == 0)", fixed = TRUE)
  expect_match(
    row,
    "(BANDWIDTH_den_extern == BW_ADAP_NN)",
    fixed = TRUE
  )
  expect_match(
    row,
    "(void)np_conditional_xrow_reciprocal_cache_try(ctx);",
    fixed = TRUE
  )
  expect_match(
    row,
    "ctx->reciprocal_cache->workspace.reciprocal_storage",
    fixed = TRUE
  )
})

test_that("admitted adaptive Gaussian CVML matches the dense signed-WLS oracle", {
  skip_if_not_installed("np")
  suppressPackageStartupMessages(library(np))
  old <- options(
    np.messages = FALSE,
    np.tree = FALSE,
    np.largeh = FALSE,
    np.macMseries.accelerate = TRUE
  )
  on.exit(options(old), add = TRUE)
  on.exit(
    tryCatch(
      .Call("C_np_shadow_reset_state", PACKAGE = "np"),
      error = function(e) NULL
    ),
    add = TRUE
  )

  set.seed(2026072903L)
  n <- 48L
  xdat <- data.frame(
    x1 = runif(n, -0.82, 0.91),
    x2 = runif(n, -0.89, 0.86)
  )
  ydat <- data.frame(
    y = sin(1.4 * xdat$x1) + 0.31 * cos(1.3 * xdat$x2) +
      rnorm(n, sd = 0.19)
  )
  bw <- npcdensbw(
    xdat = xdat,
    ydat = ydat,
    bws = c(17L, 15L, 13L),
    bandwidth.compute = FALSE,
    bwmethod = "cv.ml",
    bwtype = "adaptive_nn",
    regtype = "lp",
    basis = "glp",
    degree = c(2L, 2L),
    bernstein.basis = FALSE,
    cxkertype = "gaussian",
    cxkerorder = 8L,
    cykertype = "gaussian",
    cykerorder = 2L
  )
  empty <- matrix(numeric(), nrow = n, ncol = 0L)
  rbw <- c(
    bw$xbw[bw$ixcon],
    bw$ybw[bw$iycon],
    bw$ybw[bw$iyuno],
    bw$ybw[bw$iyord],
    bw$xbw[bw$ixuno],
    bw$xbw[bw$ixord]
  )
  result <- .Call(
    "C_np_shadow_cv_density_conditional",
    empty,
    empty,
    as.matrix(ydat),
    empty,
    empty,
    as.matrix(xdat),
    as.double(rbw),
    as.integer(2L),
    as.integer(0L),
    as.integer(0L),
    as.integer(0L),
    as.integer(3L),
    as.integer(0L),
    as.integer(0L),
    FALSE,
    as.integer(0L),
    as.integer(np:::npConditionalRegtypeCode(
      bw$regtype.engine,
      bw$degree.engine,
      bw$xncon
    )),
    as.integer(bw$degree.engine),
    as.integer(FALSE),
    as.integer(np:::npLpBasisCode(bw$basis.engine)),
    PACKAGE = "np"
  )

  expect_true(is.finite(result$proof))
  expect_true(is.finite(result$prod))
  expect_equal(result$prod, result$proof, tolerance = 5e-8)
})
