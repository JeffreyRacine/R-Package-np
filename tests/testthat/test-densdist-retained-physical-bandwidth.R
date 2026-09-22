r21_density_twins <- function(d, family, kernel, type, bounded, h) {
  args <- list(dat = d, bws = h, bwscaling = TRUE,
    bandwidth.compute = FALSE, ckertype = kernel, bwtype = type)
  if (bounded) {
    args$ckerbound <- "fixed"
    args$ckerlb <- rep(0, ncol(d))
    args$ckerub <- rep(1, ncol(d))
  }
  constructor <- if (family == "density") npudensbw else npudistbw
  scaled <- do.call(constructor, args)
  args$bws <- scaled$bandwidth$x
  args$bwscaling <- FALSE
  list(scaled = scaled, physical = do.call(constructor, args))
}

test_that("unconditional fit boundaries consume retained physical bandwidths", {
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  set.seed(612)
  d <- data.frame(x = runif(27, .02, .98), z = runif(27, .03, .97))
  e <- d[c(2, 8, 15), , drop = FALSE] + .001
  for (family in c("density", "distribution"))
    for (kernel in c("gaussian", "beta"))
      for (type in c("fixed", "generalized_nn", "adaptive_nn")) {
        fit <- if (family == "density") npudens else npudist
        twins <- r21_density_twins(d, family, kernel, type, kernel == "beta",
                                  if(type == "fixed") c(.8, .9) else c(7, 8))
        retained <- twins$scaled
        seed <- .Random.seed
        for (external in c(FALSE, TRUE)) {
          args <- list(tdat = d, se = TRUE)
          if (external) args$edat <- e
          if (type == "adaptive_nn" && kernel == "beta") {
            expect_warning(a <- do.call(fit, c(list(bws = twins$scaled), args)),
              "ANN uncertainty with finite kernel bounds is not yet implemented")
            expect_warning(b <- do.call(fit, c(list(bws = twins$physical), args)),
              "ANN uncertainty with finite kernel bounds is not yet implemented")
          } else {
            a <- do.call(fit, c(list(bws = twins$scaled), args))
            b <- do.call(fit, c(list(bws = twins$physical), args))
          }
          expect_equal(fitted(a), fitted(b), tolerance = 2e-12)
          expect_equal(se(a), se(b), tolerance = 2e-12)
          expect_equal(predict(a, newdata = e),
                       predict(b, newdata = e), tolerance = 2e-12)
        }
        expect_identical(twins$scaled, retained)
        expect_identical(.Random.seed, seed)
      }
})

test_that("fixed scaled beta estimates match independent physical contributions", {
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  d <- data.frame(x = seq(.05, .95, length.out = 24),
                  z = .5 + .4*sin(seq_len(24)))
  for (family in c("density", "distribution")) {
    fit <- if (family == "density") npudens else npudist
    b <- r21_density_twins(d, family, "beta", "fixed", TRUE, c(.8, .9))$scaled
    for (external in c(FALSE, TRUE)) {
      e <- if (external) d[c(3, 10, 19), , drop = FALSE] + .002 else d
      K <- vapply(seq_len(nrow(e)), function(i) {
        pieces <- vapply(seq_len(ncol(d)), function(j) {
          tau <- 1/b$bandwidth$x[j]^2
          if (family == "density")
            dbeta(d[[j]], 1+e[i,j]*tau, 1+(1-e[i,j])*tau) else
            pbeta(e[i,j], 1+d[[j]]*tau, 1+(1-d[[j]])*tau)
        }, numeric(nrow(d)))
        apply(pieces, 1L, prod)
      }, numeric(nrow(d)))
      args <- list(bws = b, tdat = d, se = TRUE)
      if (external) args$edat <- e
      actual <- do.call(fit, args)
      expected <- colMeans(K)
      expected.se <- sqrt(colSums(sweep(K, 2L, expected)^2)/
                            (nrow(d)*(nrow(d)-1)))
      expect_equal(fitted(actual), expected, tolerance = 2e-12)
      expect_equal(se(actual), expected.se, tolerance = 2e-12)
    }
  }
})

test_that("physical fit payloads retain categorical and bounded legacy behavior", {
  old <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old), add = TRUE)
  d <- data.frame(x = seq(.05, .95, length.out = 30),
    o = ordered(rep(0:2, 10)), u = factor(rep(letters[1:2], 15)))
  for (family in c("density", "distribution"))
    for (type in c("fixed", "generalized_nn", "adaptive_nn"))
      for (kernel in c("gaussian", "beta")) {
        dd <- if (family == "distribution") d[c("x","o")] else d
        fit <- if (family == "density") npudens else npudist
        twins <- r21_density_twins(dd, family, kernel, type, TRUE,
          c(if(type == "fixed") .6 else 7, rep(.3, ncol(dd)-1L)))
        for (tree in if(kernel == "beta") list(FALSE) else list(FALSE, TRUE, "auto")) {
          options(np.tree = tree)
          if (type == "adaptive_nn") {
            expect_warning(a <- fit(twins$scaled, tdat = dd, se = TRUE),
              "ANN uncertainty with finite kernel bounds is not yet implemented")
            expect_warning(b <- fit(twins$physical, tdat = dd, se = TRUE),
              "ANN uncertainty with finite kernel bounds is not yet implemented")
          } else {
            a <- fit(twins$scaled, tdat = dd, se = TRUE)
            b <- fit(twins$physical, tdat = dd, se = TRUE)
          }
          expect_equal(fitted(a), fitted(b), tolerance = 2e-12)
          expect_equal(se(a), se(b), tolerance = 2e-12)
        }
      }
})
