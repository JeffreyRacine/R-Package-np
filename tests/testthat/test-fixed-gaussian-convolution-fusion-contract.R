test_that("fixed ordinary-Gaussian convolution fusion is narrowly gated", {
  source_file <- test_path("..", "..", "src", "jksum.c")
  skip_if_not(file.exists(source_file), "package C sources unavailable")

  text <- paste(readLines(source_file, warn = FALSE), collapse = "\n")

  expect_match(
    text,
    "static int NP_NOINLINE np_accel_gauss_convolution_product_try",
    fixed = TRUE
  )
  expect_match(text, "(!np_accel_gauss_resolve())", fixed = TRUE)
  expect_match(text, "(BANDWIDTH_reg == BW_FIXED)", fixed = TRUE)
  expect_match(text, "(num_reg_continuous >= 2)", fixed = TRUE)
  expect_match(text, "(num_xt >= 256)", fixed = TRUE)
  expect_match(text, "(p_nvar == 0)", fixed = TRUE)
  expect_match(text, "(!int_cker_bound_extern)", fixed = TRUE)
  expect_match(
    text,
    "(operator[i] != OP_CONVOLUTION) || (KERNEL_reg[i] != 0)",
    fixed = TRUE
  )
  expect_match(
    text,
    "if(fused_gaussian_convolution_eligible)",
    fixed = TRUE
  )
  expect_match(
    text,
    "np_convol_ckernelv(KERNEL_reg[i]",
    fixed = TRUE
  )

  gate_end <- grep(
    "fused_gaussian_convolution_eligible = 1;",
    src <- readLines(source_file, warn = FALSE),
    fixed = TRUE
  )
  expect_length(gate_end, 1L)
  gate <- paste(src[seq.int(gate_end - 12L, gate_end)], collapse = "\n")
  expect_match(gate, "(BANDWIDTH_reg == BW_FIXED)", fixed = TRUE)
  expect_false(grepl("BW_GEN_NN", gate, fixed = TRUE))
  expect_false(grepl("BW_ADAP_NN", gate, fixed = TRUE))
})

test_that("fused fixed Gaussian convolution has the legacy product algebra", {
  xt <- cbind(
    c(-1.25, -0.1, 0.5, 1.75),
    c(0.8, -0.7, 1.1, 0.25),
    c(-0.2, 0.4, 1.5, -1.1)
  )
  x <- c(0.2, -0.35, 0.75)
  h <- c(0.31, 0.42, 0.27)
  hy <- c(0.22, 0.37, 0.41)
  power <- c(1L, 2L, 1L)
  normalizer <- 0.3989422803

  legacy <- rep(1, nrow(xt))
  exponent_sum <- rep(0, nrow(xt))
  scale_product <- 1
  for (j in seq_along(h)) {
    h2 <- h[j]^2 + hy[j]^2
    legacy <- legacy * (
      normalizer * h[j] * hy[j] *
        exp(-0.5 * (x[j] - xt[, j])^2 / h2) /
        (sqrt(h2) * hy[j]^power[j])
    )
    exponent_sum <- exponent_sum -
      0.5 * (x[j] - xt[, j])^2 / h2
    scale_product <- scale_product *
      normalizer * h[j] * hy[j] /
      (sqrt(h2) * hy[j]^power[j])
  }

  fused <- scale_product * exp(exponent_sum)
  expect_lt(max(abs(fused - legacy)), 2e-15)
})
