test_that("density equality freezes the union support for every contraction", {
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  x <- data.frame(u = factor(rep(c("a", "b"), c(8, 12))))
  y <- data.frame(u = factor(rep(c("a", "b", "c"), c(11, 8, 1))))
  kernel <- function(a, b, categories = 3) {
    same <- outer(as.character(a), as.character(b), "==")
    ifelse(same, .8, .2/(categories - 1))
  }
  statistic <- function(a, b) {
    A <- kernel(a, a); D <- kernel(b, b); C <- kernel(a, b)
    diag(A) <- diag(D) <- 0
    n <- length(a); m <- length(b)
    I <- sum(A)/(n*(n-1)) + sum(D)/(m*(m-1)) - 2*sum(C)/(n*m)
    V <- 2*(sum(A^2)/(n^2*(n-1)^2) + sum(D^2)/(m^2*(m-1)^2) +
              2*sum(C^2)/(n^2*m^2))
    c(In = I, Tn = I/sqrt(V))
  }
  observed <- statistic(x$u, y$u)
  set.seed(42)
  pool <- c(as.character(x$u), as.character(y$u))
  draws <- replicate(9, statistic(pool[sample.int(40, 20, TRUE)],
                                 pool[sample.int(40, 20, TRUE)]))
  before.x <- serialize(x, NULL); before.y <- serialize(y, NULL)
  bx <- npudensbw(dat = x, bws = .2, bandwidth.compute = FALSE)
  by <- npudensbw(dat = y, bws = .2, bandwidth.compute = FALSE)
  before.bw <- serialize(bx, NULL)
  results <- list(
    npdeneqtest(x, y, bw.x = .2, B = 9),
    npdeneqtest(x, y, bw.y = .2, B = 9),
    npdeneqtest(x, y, bw.x = bx, B = 9),
    npdeneqtest(x, y, bw.x = bx, bw.y = by, B = 9))
  for (result in results) {
    expect_equal(c(In = result$In, Tn = result$Tn), observed, tolerance = 2e-12)
    expect_equal(result$In.bootstrap, unname(draws["In", ]), tolerance = 2e-12)
    expect_equal(result$Tn.bootstrap, unname(draws["Tn", ]), tolerance = 2e-12)
    expect_identical(result$In.P, mean(draws["In", ] >= observed["In"]))
    expect_identical(result$Tn.P, mean(draws["Tn", ] >= observed["Tn"]))
  }
  expect_identical(serialize(x, NULL), before.x)
  expect_identical(serialize(y, NULL), before.y)
  expect_identical(serialize(bx, NULL), before.bw)
})
