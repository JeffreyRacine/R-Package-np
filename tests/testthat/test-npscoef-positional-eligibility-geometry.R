test_that("smooth coefficient eligibility uses positional evaluation columns", {
  old <- options(np.messages = FALSE); on.exit(options(old), add = TRUE)
  n <- 36L; x <- cos(seq_len(n)); z <- seq(.04, .96, length.out = n)
  y <- .4 + x * (1 + z) + sin(z * 4)
  ex <- c(-.7, .1, .6); ez <- c(.15, .4, .85)
  for (rt in c("lc", "ll", "lp")) for (type in c("fixed", "generalized_nn", "adaptive_nn")) {
    args <- list(xdat = x, zdat = z, ydat = y,
      bws = if (type == "fixed") .4 else 16, regtype = rt,
      bwtype = type, bandwidth.compute = FALSE)
    if (rt == "lp") args$degree <- 2L
    b <- do.call(npscoefbw, args)
    named.x <- setNames(data.frame(ex), b$xnames)
    named.z <- setNames(data.frame(ez), b$znames)
    oracle <- npscoef(b, exdat = named.x, ezdat = named.z, se = TRUE)
    for (shape in c("vector", "renamed", "matrix")) {
      eval.x <- switch(shape, vector = ex, renamed = data.frame(other = ex), matrix = matrix(ex))
      eval.z <- switch(shape, vector = ez, renamed = data.frame(another = ez), matrix = matrix(ez))
      prior <- serialize(list(eval.x, eval.z), NULL)
      actual <- npscoef(b, exdat = eval.x, ezdat = eval.z, se = TRUE)
      expect_identical(fitted(actual), fitted(oracle))
      expect_identical(se(actual), se(oracle))
      expect_identical(serialize(list(eval.x, eval.z), NULL), prior)
      expect_identical(predict(actual, exdat = eval.x, ezdat = eval.z), fitted(oracle))
    }
  }
})

test_that("positional eligibility preserves mixed types and omitted evaluation rows", {
  old <- options(np.messages = FALSE); on.exit(options(old), add = TRUE)
  n <- 36L
  d <- data.frame(x = cos(seq_len(n)), z = seq(.04, .96, length.out = n),
                  u = factor(rep(c("a", "b"), n / 2)))
  d$y <- 1 + d$x * d$z + as.integer(d$u)
  b <- npscoefbw(y ~ x | z + u, data = d, bws = c(.4, .25), bandwidth.compute = FALSE)
  ex <- d[c(3, 8, 14), "x", drop = FALSE]
  ez <- d[c(3, 8, 14), c("z", "u")]; ez[2, 1] <- NA_real_
  oracle <- npscoef(b, exdat = ex, ezdat = ez, se = FALSE)
  names(ez) <- c("renamed.z", "renamed.u")
  actual <- npscoef(b, exdat = ex, ezdat = ez, se = FALSE)
  expect_identical(fitted(actual), fitted(oracle))
  expect_identical(actual$nobs, oracle$nobs)
})
