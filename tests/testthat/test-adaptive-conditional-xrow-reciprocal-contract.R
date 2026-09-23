adaptive_reciprocal_radius <- function(values, k) {
  # These distinct, interior-count fixtures need no tie or extended-NN policy.
  stopifnot(!anyDuplicated(values), k == as.integer(k),
            k >= 1L, k < length(values))
  vapply(seq_along(values), function(i) {
    sort(abs(values[-i] - values[i]))[k]
  }, numeric(1L))
}

adaptive_reciprocal_gaussian8 <- function(u) {
  dnorm(u) * (105 - 105*u^2 + 21*u^4 - u^6)/48
}

adaptive_signed_wls_rows <- function(xdat, bw) {
  n <- nrow(xdat)
  stopifnot(ncol(xdat) == 2L, all(bw$degree.engine == 2L),
            identical(bw$basis.engine, "glp"), !bw$bernstein.basis.engine,
            bw$cxkerorder == 8L, bw$cykerorder == 2L)
  basis <- with(xdat, cbind(1, x1, x1^2, x2, x1*x2, x2^2))
  delete_one <- matrix(0, nrow = n, ncol = n)

  for (held_out in seq_len(n)) {
    donor <- setdiff(seq_len(n), held_out)
    # Numeric npksum dispatch returns raw exported weights: its constructor
    # dots are not the default method's private divided-row API. Build the
    # independently normalized Gaussian8 product, including every donor h.
    weight <- Reduce(`*`, lapply(seq_len(ncol(xdat)), function(j) {
      h <- adaptive_reciprocal_radius(xdat[donor, j], bw$xbw[j])
      adaptive_reciprocal_gaussian8(
        (xdat[donor, j] - xdat[held_out, j])/h)/h
    }))
    donor_basis <- basis[donor, , drop = FALSE]
    coefficient <- solve(
      crossprod(donor_basis, donor_basis * weight), basis[held_out, ]
    )
    delete_one[donor, held_out] <-
      as.numeric(weight * (donor_basis %*% coefficient))
  }
  delete_one
}

adaptive_cvml_signed_wls_oracle <- function(xdat, ydat, bw) {
  n <- nrow(xdat)
  delete_one <- adaptive_signed_wls_rows(xdat, bw)
  fit <- vapply(seq_len(n), function(held_out) {
    donor <- setdiff(seq_len(n), held_out)
    h <- adaptive_reciprocal_radius(ydat[donor, 1L], bw$ybw[1L])
    ykernel <- dnorm((ydat[donor, 1L] - ydat[held_out, 1L])/h)/h
    sum(delete_one[donor, held_out] * ykernel)
  }, numeric(1L))

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
  objective <- 0

  for (held_out in seq_len(n)) {
    donor <- setdiff(seq_len(n), held_out)
    h <- adaptive_reciprocal_radius(ydat[donor, 1L], bw$ybw[1L])
    yintegral <- vapply(ydat[[1L]], function(at) {
      pnorm((at - ydat[donor, 1L])/h)
    }, numeric(length(donor)))
    fit <- colSums(delete_one[donor, held_out] * yintegral)
    indicator <- as.numeric(ydat[[1L]][[held_out]] <= ydat[[1L]])
    keep <- seq_len(n) != held_out
    objective <- objective + sum((indicator[keep] - fit[keep])^2)
  }
  objective / (n * (n - 1L))
}

test_that("admitted adaptive Gaussian CVML retains the signed-WLS objective oracle", {
  skip_if_not_installed("npRmpi")
  skip_if_not(.mpi_pool_active(), "requires a live MPI pool")
  suppressPackageStartupMessages(library(npRmpi))
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
  objective <- npRmpi:::.npcdensbw_eval_only(xdat, ydat, bw)$objective
  oracle <- adaptive_cvml_signed_wls_oracle(xdat, ydat, bw)

  expect_true(is.finite(objective))
  expect_equal(as.numeric(objective), oracle, tolerance = 5e-8)

  held <- 7L
  donor <- setdiff(seq_len(n), held)
  fold <- npcdensbw(xdat=xdat[donor, ], ydat=ydat[donor, , drop=FALSE],
    bws=c(bw$ybw, bw$xbw), bandwidth.compute=FALSE, bwtype="adaptive_nn",
    regtype="lp", degree=c(2L,2L), basis="glp", bernstein.basis=FALSE,
    cxkertype="gaussian", cxkerorder=8L, cykertype="gaussian", cykerorder=2L)
  fitted.fold <- fitted(npcdens(bws=fold, txdat=xdat[donor, ],
    tydat=ydat[donor, , drop=FALSE], exdat=xdat[held, ], eydat=ydat[held, , drop=FALSE]))
  h <- adaptive_reciprocal_radius(ydat[donor, 1L], bw$ybw[1L])
  expected <- sum(adaptive_signed_wls_rows(xdat, bw)[donor, held] *
    dnorm((ydat[donor, 1L] - ydat[held, 1L])/h)/h)
  expect_equal(as.numeric(fitted.fold), expected, tolerance=5e-8)
})

test_that("adaptive conditional-distribution CVLS shares the signed-WLS rows", {
  skip_if_not_installed("npRmpi")
  skip_if_not(.mpi_pool_active(), "requires a live MPI pool")
  suppressPackageStartupMessages(library(npRmpi))
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
  objective <- npRmpi:::.npcdistbw_eval_only(
    xdat, ydat, bws = bw, do.full.integral = TRUE
  )$objective
  oracle <- adaptive_cdist_signed_wls_oracle(xdat, ydat, bw)

  expect_true(is.finite(objective))
  expect_equal(as.numeric(objective), oracle, tolerance = 5e-8)

  held <- 7L
  donor <- setdiff(seq_len(n), held)
  fold <- npcdistbw(xdat=xdat[donor, ], ydat=ydat[donor, , drop=FALSE],
    bws=c(bw$ybw, bw$xbw), bandwidth.compute=FALSE, bwtype="adaptive_nn",
    regtype="lp", degree=c(2L,2L), basis="glp", bernstein.basis=FALSE,
    cxkertype="gaussian", cxkerorder=8L, cykertype="gaussian", cykerorder=2L)
  fitted.fold <- fitted(npcdist(bws=fold, txdat=xdat[donor, ],
    tydat=ydat[donor, , drop=FALSE], exdat=xdat[rep(held,n), ], eydat=ydat))
  h <- adaptive_reciprocal_radius(ydat[donor, 1L], bw$ybw[1L])
  influence <- adaptive_signed_wls_rows(xdat, bw)[donor, held]
  expected <- vapply(ydat[[1L]], function(at) {
    sum(influence * pnorm((at - ydat[donor, 1L])/h))
  }, numeric(1L))
  expect_equal(as.numeric(fitted.fold), expected, tolerance=5e-8)
})

locate_mpi_conditional_xrow_source <- function() {
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

test_that("MPI adaptive conditional reciprocals remain an isolated sidecar", {
  path <- locate_mpi_conditional_xrow_source()
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
  expect_match(source, "(np_glp_cv_cache.nterms < 4)", fixed = TRUE)
  expect_match(source, "(num_train <= 0) || (ndim < 2)", fixed = TRUE)
  expect_match(
    source,
    "reciprocal_count = ((size_t)ndim + 1)*(size_t)num_train;",
    fixed = TRUE
  )
  expect_match(
    source,
    "cache->product_reciprocal =",
    fixed = TRUE
  )
  expect_match(
    source,
    "np_accel_gauss_adaptive_higher_row_reciprocal_try(",
    fixed = TRUE
  )
  guarded_helper <- paste(
    "#if NP_ACCEL_GAUSS_COMPILED",
    "static int NP_NOINLINE",
    "np_accel_gauss_adaptive_higher_row_reciprocal_try(",
    sep = "\n"
  )
  guarded_occurrences <- gregexpr(
    guarded_helper,
    source,
    fixed = TRUE
  )[[1L]]
  expect_length(guarded_occurrences[guarded_occurrences > 0L], 2L)

  prepare_start <- regexpr(
    "static int np_conditional_xrow_ctx_prepare",
    source,
    fixed = TRUE
  )
  row_start <- regexpr(
    "static int NP_NOINLINE NP_HOT_ALIGN np_conditional_xrow_from_ctx_impl",
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
    "(void)np_conditional_xrow_reciprocal_cache_try(ctx);",
    fixed = TRUE
  )
  expect_match(
    row,
    "ctx->reciprocal_cache->storage",
    fixed = TRUE
  )
  expect_match(
    row,
    "ctx->reciprocal_cache->product_reciprocal",
    fixed = TRUE
  )
})

test_that("MPI regression retains its original division-only Gaussian helper", {
  path <- locate_mpi_conditional_xrow_source()
  skip_if(is.null(path), "package C sources unavailable in this test context")
  source <- paste(readLines(path, warn = FALSE), collapse = "\n")

  helper_start <- regexpr(
    "static int NP_NOINLINE np_accel_gauss_adaptive_higher_row_try(",
    source,
    fixed = TRUE
  )
  helper_end <- regexpr(
    "static int NP_NOINLINE np_accel_gauss_product_kind(",
    source,
    fixed = TRUE
  )
  expect_gt(helper_start, 0L)
  expect_gt(helper_end, helper_start)
  helper <- substr(source, helper_start, helper_end - 1L)

  expect_false(grepl("bandwidth_reciprocal", helper, fixed = TRUE))
  expect_match(
    helper,
    "np_accel_vdivD(bandwidth[d], 1, np_accel_gauss_tmp, 1,",
    fixed = TRUE
  )
})
