test_that("npregivderiv owns ordinary-CDF adjoint normalization", {
  src_path <- testthat::test_path("..", "..", "R", "npregivderiv.R")
  skip_if_not(file.exists(src_path), "source R files unavailable")
  src <- paste(readLines(src_path, warn = FALSE), collapse = "\n")

  expect_match(
    src,
    "npksum\\.dots <- npksum\\.dots\\[names\\(npksum\\.dots\\) != \\\"bandwidth\\.divide\\\"\\]",
    perl = TRUE
  )
  expect_match(src, "bandwidth\\.divide=TRUE", perl = TRUE)
  expect_length(gregexpr(
    "cdf\\.weighted\\.average <- cdf\\.weighted\\.average\\.apply",
    src,
    perl = TRUE
  )[[1L]], 2L)
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

  fitted_second_term <- gregexpr(
    "S\\.z\\*mean\\.predicted\\.E\\.mu\\.w",
    src,
    perl = TRUE
  )[[1L]]

  expect_length(fitted_second_term[fitted_second_term > 0L], 2L)
  expect_false(grepl("S\\.z\\*mean\\.mu", src, perl = TRUE))
  expect_false(grepl("mean\\.mu <- mean\\(mu\\)", src, perl = TRUE))
})
