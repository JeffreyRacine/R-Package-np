test_that("copula marginal bounds follow original variable identity", {
  old <- options(np.messages = FALSE); on.exit(options(old), add = TRUE)
  helper <- getFromNamespace(".npcopula_marginal_bw_args", "np")
  d <- data.frame(o = ordered(rep(1:3, 8)), x = seq(.1, .9, length.out = 24),
                  z = 12 + sin(seq_len(24)))
  lo <- c(-Inf, 0, 10); hi <- c(Inf, 1, 14)
  for (type in c("fixed", "generalized_nn", "adaptive_nn"))
    for (target in c("density", "distribution"))
      for (perm in list(1:3, c(2, 3, 1), c(3, 1, 2))) {
        dd <- d[perm]
        widths <- if (type == "fixed") c(.2, .2, .7) else c(.2, 12, 12)
        constructor <- if (target == "density") npudensbw else npudistbw
        b <- do.call(constructor, list(dat = dd, bws = widths[perm], bwtype = type,
          ckerbound = "fixed", ckerlb = lo[perm], ckerub = hi[perm], bandwidth.compute = FALSE))
        for (j in seq_along(dd)) {
          a <- helper(b, dd, j, target)
          expect_identical(unname(a$bws), unname(b$bandwidth$x[j]))
          expect_identical(a$ckerlb, if (is.numeric(dd[[j]])) b$ckerlb[j] else NULL)
          expect_identical(a$ckerub, if (is.numeric(dd[[j]])) b$ckerub[j] else NULL)
        }
      }
})

test_that("bounded mixed copula sample fits agree with literal marginals", {
  old <- options(np.messages = FALSE); on.exit(options(old), add = TRUE)
  d <- data.frame(o = ordered(rep(1:3, 8)), x = seq(.1, .9, length.out = 24))
  for (kernel in c("gaussian", "epanechnikov")) for (target in c("density", "distribution")) {
    constructor <- if (target == "density") npudensbw else npudistbw
    b <- do.call(constructor, list(dat = d, bws = c(.2, .25), ckertype = kernel,
      ckerbound = "fixed", ckerlb = 0, ckerub = 1, bandwidth.compute = FALSE))
    actual <- npcopula(b, data = d, se = TRUE)
    cx <- npudistbw(dat = d["x"], bws = .25, ckertype = kernel,
      ckerbound = "fixed", ckerlb = 0, ckerub = 1, bandwidth.compute = FALSE)
    co <- npudistbw(dat = d["o"], bws = .2, bandwidth.compute = FALSE)
    coords <- as.data.frame(actual)
    expect_equal(as.numeric(coords[["u2"]]), as.numeric(fitted(npudist(cx))), tolerance = 0)
    expect_equal(as.numeric(coords[["u1"]]), as.numeric(fitted(npudist(co))), tolerance = 0)
    expected <- if (target == "distribution") fitted(npudist(b)) else {
      dx <- npudensbw(dat = d["x"], bws = .25, ckertype = kernel,
        ckerbound = "fixed", ckerlb = 0, ckerub = 1, bandwidth.compute = FALSE)
      do <- npudensbw(dat = d["o"], bws = .2, bandwidth.compute = FALSE)
      fitted(npudens(b)) / fitted(npudens(do)) / fitted(npudens(dx))
    }
    expect_equal(as.numeric(fitted(actual)), as.numeric(expected), tolerance = 1e-13)
    reverse <- do.call(constructor, list(dat = d[2:1], bws = c(.25, .2), ckertype = kernel,
      ckerbound = "fixed", ckerlb = 0, ckerub = 1, bandwidth.compute = FALSE))
    reversed <- npcopula(reverse, data = d[2:1], se = TRUE)
    expect_equal(as.numeric(fitted(actual)), as.numeric(fitted(reversed)), tolerance = 1e-13)
    expect_equal(as.numeric(se(actual)), as.numeric(se(reversed)), tolerance = 1e-13)
  }
})

test_that("beta continuous margins retain their distinct original bounds", {
  old <- options(np.messages = FALSE); on.exit(options(old), add = TRUE)
  helper <- getFromNamespace(".npcopula_marginal_bw", "np")
  d <- data.frame(o = ordered(rep(1:3, 8)), x = seq(.1, .9, length.out = 24),
                  z = 12 + sin(seq_len(24)))
  for (target in c("density", "distribution")) {
    constructor <- if (target == "density") npudensbw else npudistbw
    b <- do.call(constructor, list(dat = d, bws = c(.2, .15, .4),
      ckertype = "beta", ckerbound = "fixed", ckerlb = c(-Inf, 0, 10),
      ckerub = c(Inf, 1, 14), bandwidth.compute = FALSE))
    for (j in 2:3) {
      marginal <- helper(b, d, j, target)
      expect_identical(unname(marginal$ckerlb), c(0, 10)[j - 1L])
      expect_identical(unname(marginal$ckerub), c(1, 14)[j - 1L])
      fit <- if (target == "density") npudens(marginal) else npudist(marginal)
      expect_true(all(is.finite(fitted(fit))))
    }
  }
})
