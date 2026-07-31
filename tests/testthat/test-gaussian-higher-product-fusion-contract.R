library(npRmpi)

test_that("MPI higher-order Gaussian product fusion is narrow and rank-local", {
  src <- readLines(npRmpi_test_source_path("src", "jksum.c"), warn = FALSE)
  text <- paste(src, collapse = "\n")

  expect_match(
    text, "static int NP_NOINLINE np_accel_gauss_higher_product_try",
    fixed = TRUE
  )
  expect_match(
    text,
    "(kernel[d] < 1) || (kernel[d] > 3)",
    fixed = TRUE
  )
  expect_match(text, "np_accel_gauss_resolve()", fixed = TRUE)
  expect_match(
    text,
    "np_accel_gauss_scratch_ensure(num_xt)",
    fixed = TRUE
  )
  expect_match(
    text,
    "static int NP_NOINLINE np_accel_gauss_product_kind",
    fixed = TRUE
  )
  expect_match(
    text,
    "else if((fused_gaussian_product_eligible == 2) && (!any_cont_largeh))",
    fixed = TRUE
  )
  expect_match(
    text,
    paste0(
      "np_accel_gauss_higher_product_try(\n",
      "            KERNEL_reg_np, xtc, xc, m, num_reg_continuous"
    ),
    fixed = TRUE
  )
  expect_match(
    text,
    "np_ckernelv(KERNEL_reg_np[i]",
    fixed = TRUE
  )
})

test_that("MPI higher-order numerical formulas require an active pool", {
  skip_if(
    identical(environmentName(environment(npksum)), "npRmpi"),
    paste(
      "independent order-4/6/8 fixed, generalized-NN, and adaptive-NN",
      "oracles are covered by active one- and multi-slave sentinels"
    )
  )
})
