test_that("npregivderiv owns ordinary-CDF adjoint normalization", {
  src_path <- testthat::test_path("..", "..", "R", "npregivderiv.R")
  skip_if_not(file.exists(src_path), "source R files unavailable")
  src <- paste(readLines(src_path, warn = FALSE), collapse = "\n")

  expect_match(src, "\\.np_iv_deriv_adjoint_dots\\(list\\(\\.\\.\\.\\)\\)",
               perl = TRUE)
  expect_match(src, "bandwidth\\.divide=TRUE", perl = TRUE)
  expect_match(src, "ukertype=\"liracine\"", fixed = TRUE)
  expect_match(src, "okertype=\"liracine\"", fixed = TRUE)
  expect_match(src,
               "cdf\\.average <- cdf\\.weighted\\.average\\.apply\\(rhs, evaluation\\)",
               perl = TRUE)
  expect_match(src, "adjoint\\.apply\\(predicted\\.E\\.mu\\.w", perl = TRUE)
})

test_that("npregivderiv keeps the normal-reference operator route private", {
  src_path <- testthat::test_path("..", "..", "R", "npregivderiv.R")
  skip_if_not(file.exists(src_path), "source R files unavailable")
  src <- paste(readLines(src_path, warn = FALSE), collapse = "\n")

  expect_match(src,
               "npudensbw\\(dat=z, bwmethod=\"normal-reference\"\\)",
               perl = TRUE)
  expect_match(src, "npudens\\(tdat=z, bws=bw\\$bw\\)", perl = TRUE)
  expect_match(src, "npudist\\(tdat=z, bws=bw\\$bw\\)", perl = TRUE)
  expect_match(src, "bws=bw\\$bw,\\s*bandwidth.divide=TRUE", perl = TRUE)
  expect_false(grepl(
    "npudensbw\\(dat=z, bwmethod=\"normal-reference\",\\s*\\.\\.\\.",
    src, perl = TRUE
  ))
})

test_that("npregivderiv adjoint dots remove only operator-owned names", {
  filter.dots <- getFromNamespace(".np_iv_deriv_adjoint_dots", "npRmpi")
  dots <- structure(
    list(1, "user-u", 2, "user-o", 3, FALSE, 4, 5),
    names = c("", "ukertype", "other", "okertype", NA_character_,
              "bandwidth.divide", "other", "")
  )

  observed <- filter.dots(dots)
  expected <- dots[c(1L, 3L, 5L, 7L, 8L)]

  expect_identical(observed, expected)
  expect_identical(filter.dots(unname(dots)), unname(dots))
})

test_that("npregivderiv monotonicity guard uses only computed norms", {
  src_path <- testthat::test_path("..", "..", "R", "npregivderiv.R")
  skip_if_not(file.exists(src_path), "source R files unavailable")
  src <- paste(readLines(src_path, warn = FALSE), collapse = "\n")

  expect_match(
    src,
    "is\\.monotone\\.increasing\\(norm\\.stop\\[seq_len\\(N\\)\\]\\)",
    perl = TRUE
  )
  expect_false(grepl(
    "!is\\.monotone\\.increasing\\(norm\\.stop\\)\\s*\\{",
    src,
    perl = TRUE
  ))
})

test_that("npregivderiv centers both adjoint terms on the fitted residual", {
  src_path <- testthat::test_path("..", "..", "R", "npregivderiv.R")
  skip_if_not(file.exists(src_path), "source R files unavailable")
  src <- paste(readLines(src_path, warn = FALSE), collapse = "\n")

  expect_match(src, "rhs\\.mean <- mean\\(rhs\\)", perl = TRUE)
  expect_match(src, "survivor\\.average <- rhs\\.mean - cdf\\.average", perl = TRUE)
  expect_match(src, "survivor\\*rhs\\.mean", perl = TRUE)
  expect_match(src, "adjoint\\.apply\\(predicted\\.E\\.mu\\.w", perl = TRUE)
  expect_false(grepl("mean\\.mu <- mean\\(mu\\)", src, perl = TRUE))
})
