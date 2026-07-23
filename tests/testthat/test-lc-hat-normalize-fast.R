test_that("native LC hat normalization is byte-exact", {
  cases <- list(
    matrix(2.0, nrow = 1L, ncol = 1L),
    matrix(seq_len(15L) / 7, nrow = 3L, ncol = 5L),
    matrix(sin(seq_len(665L)), nrow = 35L, ncol = 19L)
  )

  for (kw in cases) {
    denominator <- colSums(kw)
    denominator[denominator == 0.0] <- .Machine$double.xmin
    reference <- sweep(t(kw), 1L, denominator, "/", check.margin = FALSE)
    candidate <- np:::.np_lc_hat_normalize(kw, denominator)
    expect_identical(candidate, reference)
  }
})

test_that("native LC hat normalization retains caller floor policy", {
  kw <- matrix(c(0, 0, 1, 2, -1, 1), nrow = 2L)
  denominator <- colSums(kw)
  denominator[denominator == 0.0] <- .Machine$double.xmin

  expect_identical(
    np:::.np_lc_hat_normalize(kw, denominator),
    sweep(t(kw), 1L, denominator, "/", check.margin = FALSE)
  )
})
