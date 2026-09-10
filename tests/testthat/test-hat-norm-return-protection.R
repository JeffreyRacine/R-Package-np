test_that("both native LP hat return shapes retain values and recover after errors", {
  x <- seq(-1, 1, length.out = 7L)
  design <- cbind(1, x, x^2)
  evaluation <- rbind(c(1, .2, .04), c(0, 1, .4))
  ordinary <- cbind(exp(-x^2), exp(-(x - .2)^2))
  ridged <- matrix(c(1, rep(0, 6L)), 7L, 2L)
  call <- function(symbol, weights) .Call(symbol, weights, design, evaluation,
                                         PACKAGE = "np")
  for (weights in list(ordinary, ridged)) {
    matrix.only <- call("C_np_reghat_lp_matrix_fast", weights)
    with.norm <- call("C_np_reghat_lp_matrix_norm", weights)
    expect_identical(names(with.norm), c("hat", "norm"))
    expect_identical(with.norm$hat, matrix.only)
    expect_identical(dim(with.norm$norm), c(2L, 3L))
    expect_equal(with.norm$norm[, 1L] * sqrt(with.norm$norm[, 2L]),
                 sqrt(rowSums(matrix.only^2)), tolerance = 1e-14)
    expect_identical(with.norm$norm[, 3L], c(0, 0))
  }
  invalid <- ordinary
  invalid[1L, 1L] <- NaN
  for (symbol in c("C_np_reghat_lp_matrix_fast", "C_np_reghat_lp_matrix_norm")) {
    expected <- call(symbol, ordinary)
    expect_error(call(symbol, invalid), "non-finite system", fixed = TRUE)
    expect_identical(call(symbol, ordinary), expected)
  }
})
