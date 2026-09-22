test_that("regression fits use retained physical smoothing parameters", {
  old <- options(np.messages = FALSE, np.tree = "auto")
  on.exit(options(old), add = TRUE)
  set.seed(707)
  n <- 32L
  x <- data.frame(x = rnorm(n), u = factor(rep(c("a", "b"), n/2)),
                  z = runif(n))
  y <- sin(x$x) + x$z + .2 * (x$u == "b") + rnorm(n, sd = .2)
  for (type in c("fixed", "generalized_nn", "adaptive_nn"))
    for (reg in c("lc", "ll", "lp")) {
      args <- list(xdat = x, ydat = y,
        bws = c(if(type == "fixed") 1.2 else 18, .2,
                if(type == "fixed") .9 else 20),
        bandwidth.compute = FALSE, bwscaling = TRUE, bwtype = type,
        regtype = reg, degree = if(reg == "lp") c(2L, 2L) else NULL)
      scaled <- do.call(npregbw, args)
      args$bws <- scaled$bandwidth$x
      args$bwscaling <- FALSE
      physical <- do.call(npregbw, args)
      saved <- serialize(scaled, NULL)
      for (tree in list(FALSE, TRUE, "auto")) {
        options(np.tree = tree)
        for (external in c(FALSE, TRUE)) {
          inputs <- if (external) list(exdat = x[c(3, 8, 17, 25), ]) else list()
          a <- do.call(npreg, c(list(bws = scaled, txdat = x, tydat = y,
                                    gradients = TRUE, se = TRUE), inputs))
          b <- do.call(npreg, c(list(bws = physical, txdat = x, tydat = y,
                                    gradients = TRUE, se = TRUE), inputs))
          for (field in c("mean", "merr", "grad", "gerr"))
            expect_identical(a[[field]], b[[field]])
        }
      }
      expect_identical(serialize(scaled, NULL), saved)
    }
})

test_that("scaled regression agrees with an independent Gaussian mean", {
  old <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old), add = TRUE)
  x <- data.frame(x = c(-1, -.7, -.2, .15, .3, .8, 1.4))
  y <- c(1, 0, -2, 3, 2, -1, 4)
  e <- data.frame(x = c(-.3, .1, .6))
  bw <- npregbw(xdat = x, ydat = y, bws = 1, bwscaling = TRUE,
                 bandwidth.compute = FALSE)
  k <- dnorm(outer(x$x, e$x, "-") / bw$bandwidth$x)
  oracle <- colSums(k * y) / colSums(k)
  expect_equal(fitted(npreg(bw, txdat = x, tydat = y, exdat = e)),
                oracle, tolerance = 1e-14)
})
