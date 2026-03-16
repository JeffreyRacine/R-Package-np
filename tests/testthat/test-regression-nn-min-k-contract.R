test_that("nonfixed regression NN bandwidths enforce k >= 2 across lc ll and lp", {
  set.seed(42)
  n <- 80
  x <- sort(rnorm(n))
  y <- x^2 + rnorm(n, sd = 0.15)
  tx <- data.frame(x = x)
  ex <- data.frame(x = seq(min(x), max(x), length.out = 19))

  expect_error(
    npregbw(xdat = tx, ydat = y, regtype = "lc", bwtype = "generalized_nn",
            bandwidth.compute = FALSE, bws = 1),
    "nearest-neighbor bandwidth must be at least 2|nearest-neighbor bandwidth must be in \\[2,"
  )

  expect_error(
    npregbw(xdat = tx, ydat = y, regtype = "ll", bwtype = "generalized_nn",
            bandwidth.compute = FALSE, bws = 1),
    "nearest-neighbor bandwidth must be at least 2|nearest-neighbor bandwidth must be in \\[2,"
  )

  expect_error(
    npregbw(xdat = tx, ydat = y, regtype = "lp", degree = 2L, bwtype = "generalized_nn",
            bandwidth.compute = FALSE, bws = 1),
    "nearest-neighbor bandwidth must be at least 2|nearest-neighbor bandwidth must be in \\[2,"
  )

  bw.ll <- npregbw(xdat = tx, ydat = y, regtype = "ll", bwtype = "generalized_nn",
                   bandwidth.compute = FALSE, bws = 2)
  fit.ll <- npreg(bws = bw.ll, exdat = ex)
  hat.ll <- npreghat(bws = bw.ll, txdat = tx, exdat = ex)
  apply.ll <- npreghat(bws = bw.ll, txdat = tx, exdat = ex, y = y, output = "apply")
  expect_equal(drop(hat.ll %*% y), drop(fit.ll$mean), tolerance = 1e-10)
  expect_equal(drop(apply.ll), drop(fit.ll$mean), tolerance = 1e-10)

  bw.auto <- npregbw(xdat = tx, ydat = y, regtype = "ll", bwtype = "generalized_nn", nmulti = 1)
  expect_gte(as.integer(bw.auto$bw[1]), 2L)

  bw.bad <- bw.ll
  bw.bad$bw[1] <- 1
  bw.bad$bandwidth$x[1] <- 1
  expect_error(npreg(bws = bw.bad), "nearest-neighbor bandwidth must be in \\[2,")
  expect_error(npreghat(bws = bw.bad, txdat = tx), "nearest-neighbor bandwidth must be in \\[2,")
})

test_that("semiparametric NN selectors honor the same floor on common-use routes", {
  set.seed(9202)
  n <- 60L
  z <- sort(runif(n))
  x <- 0.7 * z + rnorm(n, sd = 0.18)
  y <- 1 + 1.1 * x + sin(2 * pi * z) + rnorm(n, sd = 0.05)
  tx <- data.frame(x = x)
  tz <- data.frame(z = z)
  ex <- data.frame(x = seq(quantile(x, 0.2), quantile(x, 0.8), length.out = 25L))
  ez <- data.frame(z = seq(0.1, 0.9, length.out = 25L))

  bw.pl <- npplregbw(xdat = tx, zdat = tz, ydat = y,
                     regtype = "ll", bwtype = "generalized_nn", nmulti = 1)
  expect_gte(min(vapply(bw.pl$bw, function(one) min(as.integer(one$bw[one$icon])), integer(1))), 2L)
  fit.pl <- npplreg(bws = bw.pl, txdat = tx, tzdat = tz, tydat = y, exdat = ex, ezdat = ez)$mean
  apply.pl <- npplreghat(bws = bw.pl, txdat = tx, tzdat = tz, exdat = ex, ezdat = ez, y = y, output = "apply")
  matrix.pl <- drop(npplreghat(bws = bw.pl, txdat = tx, tzdat = tz, exdat = ex, ezdat = ez, output = "matrix") %*% y)
  expect_equal(drop(apply.pl), drop(fit.pl), tolerance = 1e-10)
  expect_equal(drop(matrix.pl), drop(fit.pl), tolerance = 1e-10)

  set.seed(9203)
  n <- 80L
  x1 <- runif(n, -1, 1)
  x2 <- rnorm(n)
  y <- sin(2 * (0.7 * x1 - 0.3 * x2)) + 0.2 * x1 * x2 + rnorm(n, sd = 0.05)
  tx <- data.frame(x1 = x1, x2 = x2)
  ex <- data.frame(x1 = seq(-0.75, 0.75, length.out = 20L),
                   x2 = seq(-0.8, 0.8, length.out = 20L))

  bw.si <- npindexbw(xdat = tx, ydat = y, regtype = "ll",
                     bwtype = "generalized_nn", nmulti = 1)
  expect_gte(as.integer(bw.si$bw[1]), 2L)
  fit.si <- npindex(bws = bw.si, txdat = tx, tydat = y, exdat = ex,
                    gradients = FALSE, errors = FALSE)$mean
  apply.si <- npindexhat(bws = bw.si, txdat = tx, exdat = ex, y = y, output = "apply", s = 0L)
  matrix.si <- drop(npindexhat(bws = bw.si, txdat = tx, exdat = ex, output = "matrix", s = 0L) %*% y)
  expect_equal(drop(apply.si), drop(fit.si), tolerance = 1e-10)
  expect_equal(drop(matrix.si), drop(fit.si), tolerance = 1e-10)

  set.seed(9204)
  n <- 80L
  x <- runif(n, -1, 1)
  z <- rnorm(n)
  y <- x^2 + rnorm(n, sd = 0.25 * sd(x))
  tx <- data.frame(x = x)
  tz <- data.frame(z = z)

  bw.sc <- npscoefbw(xdat = tx, zdat = tz, ydat = y, nmulti = 1, bwtype = "adaptive_nn")
  expect_gte(min(as.integer(bw.sc$bw)), 2L)
  fit.sc <- npscoef(bws = bw.sc, txdat = tx, tzdat = tz, tydat = y)$mean
  matrix.sc <- drop(npscoefhat(bws = bw.sc, txdat = tx, tzdat = tz) %*% y)
  expect_equal(drop(matrix.sc), drop(fit.sc), tolerance = 1e-10)
})
