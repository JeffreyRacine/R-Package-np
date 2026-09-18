test_that("exact wild helper never substitutes LC gradients for LP coefficients", {
  exact <- getFromNamespace(".np_wild_boot_from_regression_exact", "npRmpi")
  ns <- environment(exact)
  visits <- new.env(parent=emptyenv()); visits$n <- 0L
  trace(".np_wild_boot_from_regression_operator", where=ns,
    tracer=substitute(assign("n", get("n",envir=COUNTER)+1L,envir=COUNTER),
                      list(COUNTER=visits)), print=FALSE)
  on.exit(untrace(".np_wild_boot_from_regression_operator",where=ns),add=TRUE)
  set.seed(17122)
  x <- data.frame(x = runif(40), z = runif(40))
  y <- sin(4 * x$x) + x$z^2 + rnorm(40, sd = .1)
  ex <- x[c(2, 7, 18, 25), ]
  for (rt in c("lc", "ll", "lp")) {
    visits$n <- 0L
    b <- npregbw(xdat = x, ydat = y, bws = c(.3, .3), bandwidth.compute = FALSE,
      regtype = rt, degree = if (rt == "lp") c(2L, 2L) else NULL)
    mu <- fitted(npreg(b))
    set.seed(2217)
    a <- exact(x, ex, b, y, B = 5L, wild = "rademacher",
      fit.mean.train = mu, gradients = TRUE, slice.index = 1L)
    rng <- .Random.seed
    expect_identical(visits$n > 0L, rt == "lc")
    set.seed(2217)
    draws <- matrix(ifelse(runif(40L * 5L) <= .5, -1, 1), 40L)
    expect_identical(.Random.seed, rng)
    oracle <- t(vapply(seq_len(5L), function(j)
      gradients(npreg(b, txdat = x, tydat = mu + (y - mu) * draws[, j],
        exdat = ex, gradients = TRUE))[, 1L], numeric(4L)))
    expect_equal(unname(a$t), unname(oracle), tolerance = 1e-10)
  }
})

test_that("wild operator blocks use the same complete evaluation rows", {
  exact <- getFromNamespace(".np_wild_boot_from_regression_exact", "npRmpi")
  old <- options(np.plot.wild.apply.operator.threshold.bytes = Inf)
  on.exit(options(old))
  set.seed(17121)
  x <- data.frame(x = runif(30))
  y <- x$x + rnorm(30)
  ex <- data.frame(x = c(.2, NA, .4, .6))
  b <- npregbw(xdat = x, ydat = y, bws = .3, bandwidth.compute = FALSE)
  set.seed(121)
  a <- exact(x, ex, b, y, B = 5L, gradients = TRUE)
  set.seed(121)
  z <- exact(x, na.omit(ex), b, y, B = 5L, gradients = TRUE)
  expect_identical(a, z)
  options(np.plot.wild.apply.operator.threshold.bytes = 1)
  set.seed(121)
  c <- exact(x, ex, b, y, B = 5L, gradients = TRUE)
  expect_equal(c, z, tolerance = 1e-12)
})
