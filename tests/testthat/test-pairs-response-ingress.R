test_that("pairs validates every regression outcome before fitting", {
  for (invalid in list(factor(c("a", "b")), ordered(c("a", "b")), c(TRUE, FALSE))) {
    d <- data.frame(x = 1:2, y = invalid)
    set.seed(19)
    seed <- .Random.seed
    expect_error(np.pairs(c("x", "y"), d), "must be a numeric ", fixed = TRUE)
    expect_identical(.Random.seed, seed)
  }
})

test_that("pairs sends explicit bandwidths to each panel constructor once", {
  d <- data.frame(x = seq(-1, 1, length.out = 16))
  d$y <- sin(2 * d$x)
  p <- np.pairs(c("x", "y"), d, bws = .25,
                bandwidth.compute = FALSE)
  for (i in seq_len(nrow(p$pair_names))) {
    x <- d[[as.character(p$pair_names[i, 1])]]
    y <- d[[as.character(p$pair_names[i, 2])]]
    k <- outer(x, x, function(a, b) dnorm((a - b) / .25))
    a <- p$pair_kerns[[i]]
    if (p$pair_names[i, 1] == p$pair_names[i, 2])
      expect_equal(a$dens, rowMeans(k) / .25, tolerance = 1e-12)
    else
      expect_equal(a$mean, as.vector((k / rowSums(k)) %*% y), tolerance = 1e-12)
  }
})
