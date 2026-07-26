test_that("fixed ordinary-Gaussian fusion remains narrowly gated", {
  src <- readLines(test_path("..", "..", "src", "jksum.c"), warn = FALSE)
  text <- paste(src, collapse = "\n")

  expect_match(text, "static int np_accel_gauss_product_try", fixed = TRUE)
  expect_match(text, "(BANDWIDTH_reg == BW_FIXED)", fixed = TRUE)
  expect_match(text, "(num_reg_continuous >= 2)", fixed = TRUE)
  expect_match(text, "(num_xt >= 256)", fixed = TRUE)
  expect_match(text, "(pxl == NULL)", fixed = TRUE)
  expect_match(text, "(!cont_largeh_any_fixed)", fixed = TRUE)
  expect_match(text, "if(KERNEL_reg_np[i] != 0)", fixed = TRUE)
  expect_match(
    text,
    "if(fused_gaussian_product_eligible)",
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
