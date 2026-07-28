library(np)

locate_alllarge_fit_projection_source <- function() {
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

test_that("all-large fit BLAS route is narrow, bounded, and fallback-complete", {
  src_file <- locate_alllarge_fit_projection_source()
  skip_if(is.null(src_file), "source file src/jksum.c unavailable")
  lines <- readLines(src_file, warn = FALSE)
  start <- grep(
    "^int kernel_estimate_regression_categorical_tree_np\\(",
    lines
  )
  stop <- grep("^  if\\(estimation_shortcut_done\\)", lines)
  expect_length(start, 1L)
  stop <- stop[stop > start][1L]
  expect_true(is.finite(stop))
  body <- paste(lines[start:stop], collapse = "\n")
  compact <- gsub("[[:space:]]+", " ", body)

  expect_match(
    paste(lines, collapse = "\n"),
    "#define NP_REG_ALLLARGE_FIT_BLOCK_MAX_ROWS 8192",
    fixed = TRUE
  )
  expect_match(compact, "basis_is_contiguous = (!do_grad) &&", fixed = TRUE)
  expect_match(
    compact,
    "np_reg_alllarge_blas_profitable(num_obs_train, glp_nterms, num_obs_train)",
    fixed = TRUE
  )
  expect_match(
    compact,
    "np_reg_alllarge_blas_profitable(num_obs_eval, glp_nterms, num_obs_eval)",
    fixed = TRUE
  )
  expect_match(
    compact,
    "basis = basis_is_contiguous ? alloc_tmatd(num_obs_train, glp_nterms) : alloc_matd(num_obs_train, glp_nterms);",
    fixed = TRUE
  )
  expect_match(body, "np_blas_gram_int(", fixed = TRUE)
  expect_match(body, "np_blas_dgemv_t_int(", fixed = TRUE)
  expect_match(body, "np_blas_dgemv_n_int(", fixed = TRUE)
  expect_match(body, "np_blas_project_inverse_block_int(", fixed = TRUE)
  expect_match(
    compact,
    "use_fit_projection_blas = basis_is_contiguous;",
    fixed = TRUE
  )
  expect_match(
    compact,
    "} else { for(i = 0; i < glp_nterms; i++){ inverse_workspace.rhs[i] = 0.0;",
    fixed = TRUE
  )
  expect_match(
    compact,
    "if(basis_is_contiguous) free_tmat(basis); else free_mat(basis, glp_nterms);",
    fixed = TRUE
  )
  expect_false(grepl("num_obs_eval\\s*\\*\\s*num_obs_eval", body))
  expect_false(grepl("num_obs_train\\s*\\*\\s*num_obs_train", body))
})

test_that("all-large fit BLAS agrees with scalar fallback and dense WLS", {
  old_options <- options(
    np.messages = FALSE,
    np.tree = FALSE,
    np.largeh = TRUE
  )
  on.exit(options(old_options), add = TRUE)

  set.seed(2026071334L)
  n <- 2048L
  tx <- data.frame(x1 = runif(n, -0.8, 0.8), x2 = runif(n, -0.8, 0.8))
  y <- sin(2 * tx$x1) + tx$x2^2 + rnorm(n, sd = 0.2)

  make_fit <- function(degree, accelerate) {
    bw <- npregbw(
      xdat = tx,
      ydat = y,
      bws = c(1e8, 1e8),
      bandwidth.compute = FALSE,
      bwmethod = "cv.ls",
      bwtype = "fixed",
      regtype = "lp",
      basis = "glp",
      degree = degree,
      bernstein.basis = FALSE,
      ckertype = "epanechnikov",
      ckerorder = 2L
    )
    options(np.macMseries.accelerate = accelerate)
    npreg(bws = bw, txdat = tx, tydat = y, exdat = tx, gradients = FALSE)
  }

  active_scalar <- make_fit(c(2L, 2L), FALSE)
  active_blas <- make_fit(c(2L, 2L), TRUE)
  expect_equal(fitted(active_blas), fitted(active_scalar), tolerance = 1e-12)
  expect_equal(active_blas$merr, active_scalar$merr, tolerance = 1e-12)

  b <- cbind(
    1,
    tx$x1,
    tx$x2,
    tx$x1^2,
    tx$x1 * tx$x2,
    tx$x2^2
  )
  gram_inverse <- solve(crossprod(b))
  beta <- gram_inverse %*% crossprod(b, y)
  sigma2hat <- mean((y - mean(y))^2)
  oracle_mean <- as.numeric(b %*% beta)
  oracle_merr <- sqrt(
    pmax(0, sigma2hat * rowSums((b %*% gram_inverse) * b))
  )
  expect_equal(as.numeric(fitted(active_blas)), oracle_mean, tolerance = 1e-12)
  expect_equal(as.numeric(active_blas$merr), oracle_merr, tolerance = 1e-12)

  inactive_scalar <- make_fit(c(1L, 1L), FALSE)
  inactive_blas <- make_fit(c(1L, 1L), TRUE)
  expect_identical(fitted(inactive_blas), fitted(inactive_scalar))
  expect_identical(inactive_blas$merr, inactive_scalar$merr)
})
