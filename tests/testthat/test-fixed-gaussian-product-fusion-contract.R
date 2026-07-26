test_that("ordinary-Gaussian fusion admits fixed and generalized NN only", {
  src <- readLines(test_path("..", "..", "src", "jksum.c"), warn = FALSE)
  text <- paste(src, collapse = "\n")

  expect_match(text, "static int np_accel_gauss_product_try", fixed = TRUE)
  expect_match(
    text,
    "((BANDWIDTH_reg == BW_FIXED) || (BANDWIDTH_reg == BW_GEN_NN))",
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
  expect_match(text, "if(KERNEL_reg_np[i] != 0)", fixed = TRUE)
  expect_match(
    text,
    "if(fused_gaussian_product_eligible && (!any_cont_largeh))",
    fixed = TRUE
  )
  expect_false(grepl(
    paste0(
      "((BANDWIDTH_reg == BW_FIXED) || ",
      "(BANDWIDTH_reg == BW_GEN_NN) || ",
      "(BANDWIDTH_reg == BW_ADAP_NN))"
    ),
    text,
    fixed = TRUE
  ))
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
