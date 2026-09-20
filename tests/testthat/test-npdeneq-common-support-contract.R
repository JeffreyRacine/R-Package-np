test_that("density equality rejects only active asymmetric ordered kernels", {
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  x <- data.frame(o = ordered(rep(1:3, c(12, 4, 5)), levels = 1:3))
  y <- data.frame(o = ordered(rep(1:3, c(3, 14, 5)), levels = 1:3))
  bw <- npudensbw(dat = x, bws = .5, bandwidth.compute = FALSE, okertype = "racineliyan")
  expect_true(all(is.finite(fitted(npudens(bws = bw, tdat = x)))))
  set.seed(816); seed <- .Random.seed
  expect_error(npdeneqtest(x, y, bw.x = bw, B = 9), "symmetric ordered kernel")
  expect_error(npdeneqtest(x, y, bw.y = bw, B = 9), "symmetric ordered kernel")
  expect_error(npdeneqtest(x, y, okertype = "r", B = 9), "symmetric ordered kernel")
  expect_identical(.Random.seed, seed)
  for (kernel in c("liracine", "wangvanryzin")) {
    bw <- npudensbw(dat = x, bws = .5, bandwidth.compute = FALSE, okertype = kernel)
    forward <- npdeneqtest(x, y, bw.x = bw, B = 9)
    reverse <- npdeneqtest(y, x, bw.x = bw, B = 9)
    expect_equal(forward$In, reverse$In, tolerance = 2e-12)
  }
  a <- data.frame(z = seq_len(20)/20)
  expect_s3_class(npdeneqtest(a, a + .1, bws = .4, bandwidth.compute = FALSE,
                             okertype = "racineliyan", B = 9), "deneqtest")
})

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
