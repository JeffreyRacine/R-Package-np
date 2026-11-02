test_that("pairs validates every regression outcome before fitting", {
  for (invalid in list(factor(c("a", "b")), ordered(c("a", "b")), c(TRUE, FALSE))) {
    d <- data.frame(x = 1:2, y = invalid)
    set.seed(19)
    seed <- .Random.seed
    expect_error(np.pairs(c("x", "y"), d), "must be a numeric ", fixed = TRUE)
    expect_identical(.Random.seed, seed)
  }
})
