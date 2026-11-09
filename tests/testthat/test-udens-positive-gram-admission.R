test_that("positive-Gram admission preserves the old successful path", {
  check <- .np_udens_positive_gram_check
  owner <- new.env(parent = environment(check))
  environment(check) <- owner
  owner$count <- 0L
  owner$.np_udens_beta_factor <- function(covariance) {
    owner$count <- owner$count + 1L
    .np_udens_beta_factor(covariance)
  }
  covariance <- matrix(c(2, .2, .2, 1), 2L)
  before <- serialize(covariance, NULL)
  expect_null(check(covariance))
  expect_identical(owner$count, 0L)
  expect_identical(serialize(covariance, NULL), before)
  # A test-local branch witness isolates strict inequality from LAPACK rounding.
  for (value in c(-1e-10 - 1e-15, -1e-10, -1e-10 + 1e-15)) {
    owner$eigen <- local({ v <- value; function(...) list(values = c(1, v)) })
    owner$count <- 0L
    expect_null(check(covariance))
    expect_identical(owner$count, as.integer(value < -1e-10))
    expect_identical(serialize(covariance, NULL), before)
  }
})

test_that("positive-Gram admission retains normalized factor refusals", {
  # Exactly rank-deficient positive Gram matrices at different physical scales.
  x <- cbind(sin(seq_len(17L)), cos(seq_len(17L)), seq_len(17L)/17)
  for (scale in c(1e-8, 1, 1e8)) {
    covariance <- tcrossprod(x * scale)
    before <- serialize(covariance, NULL)
    expect_null(.np_udens_positive_gram_check(covariance))
    expect_identical(serialize(covariance, NULL), before)
  }
  expect_error(.np_udens_positive_gram_check(matrix(c(1, 2, 2, 1), 2L)),
               "correlation is not positive semidefinite")
  expect_error(.np_udens_positive_gram_check(matrix(c(1, 2, 3, 1), 2L)),
               "symmetry")
  expect_error(.np_udens_positive_gram_check(diag(c(1, -1))),
               "invalid beta plot covariance")
  expect_error(.np_udens_positive_gram_check(diag(c(1, NaN))))
})
