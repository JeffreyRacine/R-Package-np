test_that("copula quasi-inverse grids respect the fitted marginal domain", {
  old <- options(np.messages = FALSE); on.exit(options(old), add = TRUE)
  d <- data.frame(x = seq(.1, .9, length.out = 24), z = .5 + .35 * sin(seq_len(24)))
  u <- data.frame(x = c(.2, .8), z = c(.3, .7))
  for (kernel in c("gaussian", "epanechnikov", "beta"))
    for (domain in c("fixed", "range")) for (target in c("density", "distribution")) {
      constructor <- if (target == "density") npudensbw else npudistbw
      args <- list(dat = d, bws = c(.18, .22), bandwidth.compute = FALSE,
                   ckertype = kernel, ckerbound = domain)
      if (domain == "fixed") { args$ckerlb <- 0; args$ckerub <- 1 }
      b <- do.call(constructor, args)
      actual <- npcopula(b, data = d, u = u, n.quasi.inv = 40L, er.quasi.inv = 1, se = TRUE)
      inverses <- vector("list", ncol(d))
      for (j in seq_along(d)) {
        extended <- extendrange(d[[j]], f = 1)
        extended <- c(max(extended[1], b$ckerlb[j]), min(extended[2], b$ckerub[j]))
        grid <- sort(c(seq(extended[1], extended[2], length.out = 20L),
                       quantile(d[[j]], seq(0, 1, length.out = 20L))))
        mb <- npudistbw(dat = d[j], bws = b$bandwidth$x[j], bandwidth.compute = FALSE,
          ckertype = kernel, ckerbound = "fixed", ckerlb = b$ckerlb[j], ckerub = b$ckerub[j])
        F <- fitted(npudist(mb, edat = data.frame(grid)))
        probabilities <- pmin(pmax(u[[j]], min(F)), max(F))
        inverses[[j]] <- vapply(probabilities, function(p) min(grid[F >= p]), numeric(1L))
      }
      expected.grid <- expand.grid(inverses); names(expected.grid) <- names(d)
      actual.grid <- as.data.frame(actual)[names(d)]
      expect_equal(unname(as.matrix(actual.grid)), unname(as.matrix(expected.grid)), tolerance = 0)
      expect_true(all(vapply(seq_along(d), function(j)
        all(actual.grid[[j]] >= b$ckerlb[j] & actual.grid[[j]] <= b$ckerub[j]), logical(1L))))
      expected <- if (target == "distribution") fitted(npudist(b, edat = expected.grid)) else {
        joint <- fitted(npudens(b, edat = expected.grid))
        for (j in seq_along(d)) {
          mb <- npudensbw(dat = d[j], bws = b$bandwidth$x[j], bandwidth.compute = FALSE,
            ckertype = kernel, ckerbound = "fixed", ckerlb = b$ckerlb[j], ckerub = b$ckerub[j])
          joint <- joint / fitted(npudens(mb, edat = expected.grid[j]))
        }
        joint
      }
      expect_equal(as.numeric(fitted(actual)), as.numeric(expected), tolerance = 1e-13)
      expect_identical(predict(actual, u = u, n.quasi.inv = 40L, er.quasi.inv = 1), fitted(actual))
      expect_true(all(is.finite(se(actual))))
    }
})

test_that("copula inverse extensions respect one-sided and unbounded domains", {
  old <- options(np.messages = FALSE); on.exit(options(old), add = TRUE)
  d <- data.frame(x = seq(.1, .9, length.out = 24), z = .5 + .35 * sin(seq_len(24)))
  for (bounds in list(c(-Inf, Inf), c(0, Inf), c(-Inf, 1), c(0, 1)))
    for (extension in c(0, 1)) {
      b <- npudistbw(dat = d, bws = c(.18, .22), bandwidth.compute = FALSE,
        ckerbound = "fixed", ckerlb = bounds[1L], ckerub = bounds[2L])
      actual <- npcopula(b, data = d, u = data.frame(x = .4, z = .6),
        n.quasi.inv = 40L, er.quasi.inv = extension)
      coordinates <- as.data.frame(actual)[names(d)]
      for (j in seq_along(d)) {
        extended <- extendrange(d[[j]], f = extension)
        extended <- c(max(extended[1L], bounds[1L]), min(extended[2L], bounds[2L]))
        grid <- sort(c(seq(extended[1L], extended[2L], length.out = 20L),
          quantile(d[[j]], seq(0, 1, length.out = 20L))))
        mb <- npudistbw(dat = d[j], bws = c(.18, .22)[j], bandwidth.compute = FALSE,
          ckerbound = "fixed", ckerlb = bounds[1L], ckerub = bounds[2L])
        F <- fitted(npudist(mb, edat = data.frame(grid)))
        p <- min(max(c(.4, .6)[j], min(F)), max(F))
        expect_identical(coordinates[[j]], min(grid[F >= p]))
      }
    }
})

test_that("bounded copula grids remain usable by bootstrap plots", {
  old <- options(np.messages = FALSE); on.exit(options(old), add = TRUE)
  d <- data.frame(x = seq(.1, .9, length.out = 24), z = .5 + .35 * sin(seq_len(24)))
  b <- npudensbw(dat = d, bws = c(.18, .22), bandwidth.compute = FALSE,
    ckerbound = "fixed", ckerlb = 0, ckerub = 1)
  fit <- npcopula(b, data = d, u = data.frame(x = c(.2, .8), z = c(.3, .7)),
    n.quasi.inv = 40L, er.quasi.inv = 1)
  set.seed(111)
  out <- plot(fit, output = "data", errors = "bootstrap", B = 3L, band = "pmzsd")
  expect_equal(nrow(out), 4L)
  expect_true(all(is.finite(as.matrix(out[c("center", "lower", "upper")]))))
  expect_true(all(out$lower <= out$upper))
})
