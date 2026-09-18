test_that("LC and LP0 wild gradients reuse exact operators without response refits", {
  withr::local_options(np.messages = FALSE)
  withr::local_preserve_seed()
  set.seed(180921)
  n <- 48L
  x <- data.frame(x = runif(n, -.7, .7), u = factor(rep(1:2, n/2)),
                  z = runif(n, -.7, .7), o = ordered(rep(1:3, n/3)))
  y <- sin(2*x$x) + x$z^2 + rnorm(n, sd = .2)
  ex <- x[c(2, 11, 20, 33, 41), , drop = FALSE]
  reuse <- getFromNamespace(".np_wild_boot_from_regression_exact", "npRmpi")
  for (type in c("fixed", "generalized_nn", "adaptive_nn"))
    for (kernel in c("gaussian", "epanechnikov", "uniform"))
      for (regtype in c("lc", "lp")) {
        b <- npregbw(xdat = x, ydat = y,
          bws = c(if (type == "fixed") .8 else 35, .2,
                  if (type == "fixed") .8 else 35, .25),
          bandwidth.compute = FALSE, bwtype = type, ckertype = kernel,
          regtype = regtype, degree = if (regtype == "lp") c(0L, 0L) else NULL)
        pilot <- fitted(npreg(b, txdat = x, tydat = y))
        for (wild in c("rademacher", "mammen")) {
          set.seed(819)
          u <- matrix(runif(n * 7L), n, 7L)
          draws <- if (wild == "rademacher") ifelse(u <= .5, -1, 1) else {
            a <- (1 - sqrt(5))/2
            ifelse(u <= (sqrt(5) + 1)/(2*sqrt(5)), a, 1-a)
          }
          seed <- .Random.seed
          expected <- t(vapply(seq_len(7L), function(j)
            gradients(npreg(b, txdat = x, tydat = pilot + (y-pilot)*draws[,j],
              exdat = ex, gradients = TRUE))[, 1L], numeric(nrow(ex))))
          expected.t0 <- gradients(npreg(b, txdat = x, tydat = y,
                                         exdat = ex, gradients = TRUE))[,1L]
          for (threshold in c(8*n*nrow(ex) - 1, 8*n*nrow(ex), 8*n*nrow(ex) + 1)) {
            local({
              withr::local_options(np.plot.wild.apply.operator.threshold.bytes = threshold,
                                  np.plot.wild.hat.block.bytes = 8*n*2L)
              testthat::local_mocked_bindings(.np_regression_direct = function(...)
                stop("unexpected per-response refit"), .package = "npRmpi")
              set.seed(819)
              actual <- reuse(xdat = x, exdat = ex, bws = b, ydat = y,
                fit.mean.train = pilot, B = 7L, wild = wild, gradients = TRUE, slice.index = 1L)
              expect_equal(actual$t, expected, tolerance = 1e-10)
              expect_equal(actual$t0, expected.t0, tolerance = 1e-10)
              expect_identical(.Random.seed, seed)
            })
          }
        }
      }
})

test_that("existing fixed-only wild selector decisions remain unchanged", {
  select <- getFromNamespace(".np_plot_wild_apply_operator_enabled", "npRmpi")
  budget <- getFromNamespace(".np_plot_wild_operator_exceeds_budget", "npRmpi")
  withr::local_options(np.plot.wild.apply.operator.threshold.bytes = 800)
  for (m in 9:11) {
    expect_identical(budget(10, m), m >= 10L)
    expect_identical(select(10, m, "fixed"), m >= 10L)
    expect_false(select(10, m, "generalized_nn"))
    expect_false(select(10, m, "adaptive_nn"))
  }
})
