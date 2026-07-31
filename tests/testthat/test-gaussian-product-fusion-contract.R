library(npRmpi)

test_that("ordinary-Gaussian fusion admits all scalar-bandwidth row owners", {
  src <- readLines(npRmpi_test_source_path("src", "jksum.c"), warn = FALSE)
  text <- paste(src, collapse = "\n")

  expect_match(
    text, "static int NP_NOINLINE np_accel_gauss_product_try",
    fixed = TRUE
  )
  expect_match(
    text,
    paste0(
      "((BANDWIDTH_reg == BW_FIXED) || (BANDWIDTH_reg == BW_GEN_NN) ||\n",
      "      (BANDWIDTH_reg == BW_ADAP_NN))"
    ),
    fixed = TRUE
  )
  expect_match(text, "(num_reg_continuous >= 2)", fixed = TRUE)
  expect_match(text, "(num_xt >= 256)", fixed = TRUE)
  expect_match(text, "(pxl == NULL)", fixed = TRUE)
  expect_match(
    text,
    "((BANDWIDTH_reg != BW_FIXED) || (!cont_largeh_any_fixed))",
    fixed = TRUE
  )
  expect_match(
    text, "((kind == 1) && (kernel[d] != 0))",
    fixed = TRUE
  )
  expect_match(
    text,
    "if((fused_gaussian_product_eligible == 1) && (!any_cont_largeh))",
    fixed = TRUE
  )
  expect_match(
    text,
    "const int num_xt = is_adaptive?num_obs_eval:num_obs_train;",
    fixed = TRUE
  )
  expect_match(
    text,
    paste0(
      "double * const * const xtc = is_adaptive?\n",
      "    matrix_X_continuous_eval:matrix_X_continuous_train;"
    ),
    fixed = TRUE
  )
  expect_match(
    text,
    paste0(
      "double * const * const xc = is_adaptive?\n",
      "    matrix_X_continuous_train:matrix_X_continuous_eval;"
    ),
    fixed = TRUE
  )
  expect_match(
    text,
    "const int jbw = (BANDWIDTH_reg != BW_FIXED) ? j:0;",
    fixed = TRUE
  )
  expect_match(
    text,
    "dband *= ipow(m[i][jbw], bpow[i]);",
    fixed = TRUE
  )
  expect_match(
    text,
    "if(fused_gaussian_product){",
    fixed = TRUE
  )
  expect_match(
    text,
    "np_ckernelv(KERNEL_reg_np[i]",
    fixed = TRUE
  )
})

test_that("adaptive Gaussian fusion matches observation-specific bandwidths", {
  skip_if(
    identical(environmentName(environment(npksum)), "npRmpi"),
    "npRmpi numerical oracle is covered by active-pool parity tests"
  )

  adaptive_radius <- function(train, k) {
    vapply(seq_along(train), function(i) {
      distance <- abs(train - train[[i]])
      positive <- sort(distance[distance > 0])
      duplicate_count <- sum(distance == 0) - 1L
      if (duplicate_count >= k) positive[[1L]] else
        positive[[max(1L, k - duplicate_count)]]
    }, numeric(1L))
  }

  set.seed(2026072909L)
  ntrain <- 260L
  neval <- 256L
  xtrain <- data.frame(
    x1 = runif(ntrain, -0.9, 0.9),
    x2 = runif(ntrain, -0.8, 0.85)
  )
  xeval <- data.frame(
    x1 = runif(neval, -0.87, 0.88),
    x2 = runif(neval, -0.78, 0.82)
  )
  k <- c(79L, 83L)
  h <- Map(adaptive_radius, xtrain, k)
  coefficient <- 1 / sqrt(2 * pi)
  oracle <- outer(
    xtrain$x1, xeval$x1,
    function(train, eval) {
      coefficient * exp(-0.5 * ((eval - train) / h[[1L]])^2)
    }
  ) * outer(
    xtrain$x2, xeval$x2,
    function(train, eval) {
      coefficient * exp(-0.5 * ((eval - train) / h[[2L]])^2)
    }
  )

  raw <- npksum(
    bws = k,
    txdat = xtrain,
    exdat = xeval,
    bwtype = "adaptive_nn",
    ckertype = "gaussian",
    ckerorder = 2L,
    return.kernel.weights = TRUE
  )
  divided <- npksum(
    bws = k,
    txdat = xtrain,
    exdat = xeval,
    bwtype = "adaptive_nn",
    ckertype = "gaussian",
    ckerorder = 2L,
    bandwidth.divide = TRUE,
    return.kernel.weights = TRUE
  )
  oracle_divided <- oracle /
    matrix(h[[1L]] * h[[2L]], nrow = ntrain, ncol = neval)

  expect_equal(as.matrix(raw$kw), oracle, tolerance = 2e-13)
  expect_equal(as.numeric(raw$ksum), colSums(oracle), tolerance = 2e-13)
  expect_equal(as.matrix(divided$kw), oracle_divided, tolerance = 2e-12)
  expect_equal(
    as.numeric(divided$ksum), colSums(oracle_divided), tolerance = 2e-12
  )
})
