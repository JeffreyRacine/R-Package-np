test_that("partial-linearity test is exported and both result methods are registered", {
  ns <- asNamespace("np")
  expect_identical(getExportedValue("np", "nppltest"), ns$nppltest)
  for (generic in c("print", "summary"))
    expect_identical(getS3method(generic, "pltest", envir = ns),
                     get(paste0(generic, ".pltest"), envir = ns, inherits = FALSE))
  expect_error(nppltest(NULL, B = TRUE), "B must")
  expect_error(nppltest(NULL), "object must")
})

test_that("mixed training operator uses the fitted categorical weights", {
  set.seed(481)
  n <- 30L
  x <- data.frame(x = rnorm(n), w = factor(rep(1:2, length.out = n)))
  z <- data.frame(group = factor(rep(1:2, each = n/2)), v = rnorm(n))
  y <- rnorm(n)
  responses <- cbind(y, cos(y), seq_len(n)/n)
  for (reg in c("lc", "ll")) {
    bw <- npplregbw(xdat = x, ydat = y, zdat = z, bandwidth.compute = FALSE,
                    regtype = reg, bws = rbind(c(.496, .8), c(.2, .55), c(.496, .6)))
    H <- .pltest_operator(bw, x, y, z)$H
    oracle <- apply(responses, 2L, function(v)
      as.numeric(fitted(npplreg(bws = bw, txdat = x, tydat = v, tzdat = z))))
    expect_equal(unname(H %*% responses), unname(oracle), tolerance = 1e-10)
    expect_equal(.pltest_residual_norms(H), rowSums((diag(n) - H)^2))
  }
})

test_that("Sobol prefixes use a fixed shift and tied empirical thresholds", {
  expect_identical(formals(nppltest)$directions, 2048L)
  expect_error(nppltest(NULL, directions = 0), "directions must")
  expect_error(nppltest(NULL, directions = 1.5), "directions must")
  expect_error(nppltest(NULL, directions = 2, integration = "exact"), "only")
  expect_error(.pltest_sobol_directions(1L, 1112L), "1111")

  set.seed(65)
  state <- .Random.seed
  V <- cbind(c(0, 0, 1, 2), c(1, 1, 0, 0))
  directions <- .pltest_sobol_directions(23L, 2L)
  expect_identical(directions[1:5, ], .pltest_sobol_directions(5L, 2L))
  expect_identical(state, .Random.seed)
  E <- cbind(c(2, -1, 3, -2), c(-1, 4, 1, 0))
  direct <- vapply(seq_len(nrow(directions)), function(j) {
    z <- drop(V %*% directions[j, ])
    colSums((outer(z, z, ">=") %*% E / nrow(E))^2)
  }, numeric(ncol(E)))
  expect_equal(.pltest_direction_scores(E, V, directions), rowMeans(direct), tolerance = 1e-13)
  expect_equal(.pltest_direction_scores(E * 1e150, V, directions) / 1e300,
               rowMeans(direct), tolerance = 1e-13)
  expect_equal(.pltest_direction_scores(E * 1e-150, V, directions) / 1e-300,
               rowMeans(direct), tolerance = 1e-13)
})
