test_that("conditional CDF pairs recalculate both denominators for resampled rows", {
  x <- c(-1, -.2, .6, 1.4)
  y <- c(.5, 1, -1, 2)
  grid <- data.frame(x = c(-.3, .5), y = c(.8, 1.6))
  h <- c(.4, .6)
  ids <- c(1L, 1L, 2L, 4L)
  ratio <- function(s, rows) {
    w <- dnorm(outer(x[rows], grid$x, "-") / (s*h[1]))
    num <- w * pnorm(-outer(y[rows], grid$y, "-") / (s*h[2]))
    list(w = w, num = num, value = colSums(num) / colSums(w))
  }
  a <- ratio(1, 1:4)
  b <- ratio(sqrt(2), 1:4)
  center <- 2*a$value - b$value
  plan <- list(ix = 1:2, w = a$w, jw = b$w, main = a$num, wide = b$num, center = center)
  got <- .np_conditional_pairs_block(list(value = matrix(tabulate(ids, 4), 4)), plan)
  expect_equal(as.numeric(got), 2*ratio(1, ids)$value - ratio(sqrt(2), ids)$value - center)
  plan$w[] <- 0
  expect_error(.np_conditional_pairs_block(list(value = matrix(1L, 4)), plan), "predictor mass")
})
