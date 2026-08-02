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

adaptive_xrow_radius <- function(values, k) {
  vapply(seq_along(values), function(index) {
    distance <- abs(values - values[[index]])
    duplicates <- sum(distance == 0) - 1L
    positive <- sort(distance[distance > 0])
    if (duplicates >= k) positive[[1L]] else
      positive[[max(1L, k - duplicates)]]
  }, numeric(1L))
}

adaptive_signed_wls_rows <- function(xdat, bw) {
  n <- nrow(xdat)
  basis <- np:::W.lp(
    xdat = xdat,
    degree = bw$degree.engine,
    basis = bw$basis.engine,
    bernstein.basis = bw$bernstein.basis.engine
  )
  xraw <- npksum(
    bws = bw$xbw,
    txdat = xdat,
    exdat = xdat,
    bwtype = bw$type,
    ckertype = bw$cxkertype,
    ckerorder = bw$cxkerorder,
    operator = "normal",
    return.kernel.weights = TRUE,
    bandwidth.divide = FALSE
  )$kw
  xdivisor <- Reduce(`*`, Map(
    adaptive_xrow_radius, xdat, as.list(bw$xbw)
  ))
  xraw <- xraw / matrix(xdivisor, nrow = n, ncol = n)
  delete_one <- matrix(0, nrow = n, ncol = n)

  for (evaluation in seq_len(n)) {
    weight <- xraw[, evaluation]
    coefficient <- solve(
      crossprod(basis, basis * weight), basis[evaluation, ]
    )
    full <- as.numeric(weight * (basis %*% coefficient))
    denominator <- 1 - full[[evaluation]]
    stopifnot(is.finite(denominator), denominator != 0)
    full[[evaluation]] <- 0
    delete_one[, evaluation] <- full / denominator
  }
  delete_one
}

adaptive_cvml_signed_wls_oracle <- function(xdat, ydat, bw) {
  n <- nrow(xdat)
  delete_one <- adaptive_signed_wls_rows(xdat, bw)
  ybandwidth <- adaptive_xrow_radius(ydat[[1L]], bw$ybw[[1L]])
  ydifference <- outer(
    ydat[[1L]], ydat[[1L]], function(training, evaluation) {
      evaluation - training
    }
  )
  ykernel <- dnorm(ydifference / matrix(ybandwidth, n, n)) /
    matrix(ybandwidth, n, n)
  fit <- colSums(delete_one * ykernel)

  contribution <- vapply(fit, function(value) {
    if (value > .Machine$double.xmin) return(-log(value))
    if (value < - .Machine$double.xmin)
      return(log(-value) - 2 * log(.Machine$double.xmin))
    -log(.Machine$double.xmin)
  }, numeric(1L))
  -sum(contribution)
}

adaptive_cdist_signed_wls_oracle <- function(xdat, ydat, bw) {
  n <- nrow(xdat)
  delete_one <- adaptive_signed_wls_rows(xdat, bw)
  ybandwidth <- adaptive_xrow_radius(ydat[[1L]], bw$ybw[[1L]])
  ydifference <- outer(
    ydat[[1L]], ydat[[1L]], function(training, evaluation) {
      evaluation - training
    }
  )
  yintegral <- pnorm(ydifference / matrix(ybandwidth, n, n))
  objective <- 0

  for (training in seq_len(n)) {
    fit <- colSums(delete_one[, training] * yintegral)
    indicator <- as.numeric(ydat[[1L]][[training]] <= ydat[[1L]])
    keep <- seq_len(n) != training
    objective <- objective + sum((indicator[keep] - fit[keep])^2)
  }
  objective / (n * n)
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

test_that("admitted adaptive Gaussian CVML retains the signed-WLS objective oracle", {
  skip_if_not_installed("np")
  suppressPackageStartupMessages(library(np))
  old <- options(
    np.messages = FALSE,
    np.tree = FALSE,
    np.largeh = FALSE,
    np.macMseries.accelerate = TRUE
  )
  on.exit(options(old), add = TRUE)

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
  objective <- np:::.npcdensbw_eval_only(xdat, ydat, bw)$objective
  oracle <- adaptive_cvml_signed_wls_oracle(xdat, ydat, bw)

  expect_true(is.finite(objective))
  expect_equal(as.numeric(objective), oracle, tolerance = 5e-8)
})

test_that("adaptive conditional-distribution CVLS shares the signed-WLS rows", {
  skip_if_not_installed("np")
  suppressPackageStartupMessages(library(np))
  old <- options(
    np.messages = FALSE,
    np.tree = FALSE,
    np.largeh = FALSE,
    np.macMseries.accelerate = TRUE
  )
  on.exit(options(old), add = TRUE)

  set.seed(2026080104L)
  n <- 48L
  xdat <- data.frame(
    x1 = runif(n, -0.82, 0.91),
    x2 = runif(n, -0.89, 0.86)
  )
  ydat <- data.frame(
    y = sin(1.4 * xdat$x1) + 0.31 * cos(1.3 * xdat$x2) +
      rnorm(n, sd = 0.19)
  )
  bw <- npcdistbw(
    xdat = xdat,
    ydat = ydat,
    bws = c(17L, 15L, 13L),
    bandwidth.compute = FALSE,
    bwmethod = "cv.ls",
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
  objective <- np:::.npcdistbw_eval_only(
    xdat, ydat, bws = bw, do.full.integral = TRUE
  )$objective
  oracle <- adaptive_cdist_signed_wls_oracle(xdat, ydat, bw)

  expect_true(is.finite(objective))
  expect_equal(as.numeric(objective), oracle, tolerance = 5e-8)
})
