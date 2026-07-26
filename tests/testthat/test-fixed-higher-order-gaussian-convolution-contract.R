test_that("fixed higher-order Gaussian convolution is an isolated fallback", {
  source_file <- test_path("..", "..", "src", "jksum.c")
  helper_file <- test_path(
    "..", "..", "src", "jksum_gaussian_fixed.c"
  )
  header_file <- test_path(
    "..", "..", "src", "jksum_gaussian_fixed.h"
  )
  skip_if_not(
    all(file.exists(c(source_file, helper_file, header_file))),
    "package C sources unavailable"
  )

  source <- paste(readLines(source_file, warn = FALSE), collapse = "\n")
  helper <- paste(readLines(helper_file, warn = FALSE), collapse = "\n")
  header <- paste(readLines(header_file, warn = FALSE), collapse = "\n")

  expect_match(
    source,
    "int fused_gaussian_fixed_higher_convolution_eligible = 0;",
    fixed = TRUE
  )
  expect_match(source, "(BANDWIDTH_reg == BW_FIXED)", fixed = TRUE)
  expect_match(source, "(operator[i] != OP_CONVOLUTION)", fixed = TRUE)
  expect_match(
    source,
    "((KERNEL_reg[i] != 1) && (KERNEL_reg[i] != 2))",
    fixed = TRUE
  )
  expect_match(
    source,
    "np_fixed_gaussian_convolution_product_try(",
    fixed = TRUE
  )
  expect_match(
    source,
    "np_convol_ckernelv(KERNEL_reg[i]",
    fixed = TRUE
  )
  expect_match(
    helper,
    "result[i] = xw[j]*kval/hy_power;",
    fixed = TRUE
  )
  expect_match(header, "attribute_hidden int", fixed = TRUE)

  direct <- regexpr(
    "if(fused_gaussian_convolution_kind != 0)",
    source,
    fixed = TRUE
  )
  fallback <- regexpr(
    "if(fused_gaussian_fixed_higher_convolution_eligible)",
    source,
    fixed = TRUE
  )
  expect_gt(direct, 0L)
  expect_gt(fallback, direct)
})

test_that("fixed higher-order Gaussian CVLS retains frozen objectives", {
  evaluate <- function(dat, bws, order) {
    bw <- npudensbw(
      dat = dat,
      bws = bws,
      bandwidth.compute = FALSE,
      bwmethod = "cv.ls",
      bwtype = "fixed",
      ckertype = "gaussian",
      ckerorder = order
    )
    value <- getFromNamespace("npudensbw.bandwidth", "np")(
      dat = dat,
      bws = bw,
      bandwidth.compute = TRUE,
      bwsolver = "powell",
      eval.only = TRUE,
      nmulti = 1L,
      powell.remin = FALSE
    )
    if (!is.null(value$objective))
      as.numeric(value$objective[[1L]])
    else
      -as.numeric(value$fval[[1L]])
  }

  set.seed(2026072501L + 257L)
  dat <- data.frame(
    x1 = rnorm(257L),
    x2 = rnorm(257L),
    x3 = rnorm(257L),
    x4 = rnorm(257L)
  )
  expected <- c(
    order4_q1 = -0.1600367782880025,
    order4_q2 = 0.023778818772456245,
    order6_q1 = -0.59980518555491569,
    order6_q2 = -0.16100514309017225,
    order8_q2 = -0.060404396873481092
  )
  observed <- c(
    order4_q1 = evaluate(dat["x1"], 0.31, 4L),
    order4_q2 = evaluate(dat[c("x1", "x2")], c(0.31, 0.35), 4L),
    order6_q1 = evaluate(dat["x1"], 0.31, 6L),
    order6_q2 = evaluate(dat[c("x1", "x2")], c(0.31, 0.35), 6L),
    order8_q2 = evaluate(dat[c("x1", "x2")], c(0.31, 0.35), 8L)
  )

  expect_equal(unname(observed), unname(expected), tolerance = 5e-14)
})
