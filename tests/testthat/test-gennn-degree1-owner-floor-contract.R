library(np)

test_that("generalized-nn ll and canonical lp degree-1 owners stay exact at NN floor", {
  set.seed(20260315)

  n <- 60L
  x <- rnorm(n)
  y <- x^2 + rnorm(n, sd = 0.2)
  tx <- data.frame(x = x)
  ex <- data.frame(x = seq(min(x), max(x), length.out = 31L))

  make_bw <- function(regtype) {
    args <- list(
      xdat = tx,
      ydat = y,
      regtype = regtype,
      bwtype = "generalized_nn",
      bws = 2L,
      bandwidth.compute = FALSE
    )
    if (identical(regtype, "lp")) {
      args$basis <- "glp"
      args$degree <- 1L
      args$bernstein.basis <- FALSE
    }
    do.call(npregbw, args)
  }

  bw.ll <- make_bw("ll")
  bw.lp <- make_bw("lp")

  fit.ll <- npreg(bws = bw.ll, txdat = tx, tydat = y, exdat = ex)
  fit.lp <- npreg(bws = bw.lp, txdat = tx, tydat = y, exdat = ex)

  H.ll <- npreghat(bws = bw.ll, txdat = tx, exdat = ex, output = "matrix")
  H.lp <- npreghat(bws = bw.lp, txdat = tx, exdat = ex, output = "matrix")
  hat.apply.ll <- npreghat(bws = bw.ll, txdat = tx, exdat = ex, y = y, output = "apply")
  hat.apply.lp <- npreghat(bws = bw.lp, txdat = tx, exdat = ex, y = y, output = "apply")

  expect_equal(as.vector(fit.ll$mean), as.vector(fit.lp$mean), tolerance = 1e-10)
  expect_equal(as.vector(hat.apply.ll), as.vector(hat.apply.lp), tolerance = 1e-10)
  expect_equal(as.vector(H.ll %*% y), as.vector(hat.apply.ll), tolerance = 1e-10)
  expect_equal(as.vector(H.lp %*% y), as.vector(hat.apply.lp), tolerance = 1e-10)
  expect_equal(as.vector(H.ll %*% y), as.vector(fit.ll$mean), tolerance = 1e-10)
  expect_equal(as.vector(H.lp %*% y), as.vector(fit.lp$mean), tolerance = 1e-10)
})

test_that("npplreg generalized-nn evaluation hats stay aligned at NN floor", {
  set.seed(20260315)

  n <- 80L
  z <- runif(n)
  x <- 0.7 * z + rnorm(n, sd = 0.2)
  y <- 1 + 1.1 * x + sin(2 * pi * z) + rnorm(n, sd = 0.05)

  tx <- data.frame(x = x)
  tz <- data.frame(z = z)
  ex <- data.frame(x = seq(min(x), max(x), length.out = 31L))
  ez <- data.frame(z = seq(min(z), max(z), length.out = 31L))
  bws.fixed <- matrix(c(2L, 9L), nrow = 2L)

  make_bw <- function(regtype) {
    args <- list(
      xdat = tx,
      zdat = tz,
      ydat = y,
      regtype = regtype,
      bwtype = "generalized_nn",
      bws = bws.fixed,
      bandwidth.compute = FALSE
    )
    if (identical(regtype, "lp")) {
      args$basis <- "glp"
      args$degree <- 1L
      args$bernstein.basis <- FALSE
    }
    do.call(npplregbw, args)
  }

  bw.ll <- make_bw("ll")
  bw.lp <- make_bw("lp")

  fit.ll <- npplreg(bws = bw.ll, txdat = tx, tzdat = tz, tydat = y, exdat = ex, ezdat = ez)
  fit.lp <- npplreg(bws = bw.lp, txdat = tx, tzdat = tz, tydat = y, exdat = ex, ezdat = ez)

  H.ll <- npplreghat(bws = bw.ll, txdat = tx, tzdat = tz, exdat = ex, ezdat = ez, output = "matrix")
  H.lp <- npplreghat(bws = bw.lp, txdat = tx, tzdat = tz, exdat = ex, ezdat = ez, output = "matrix")
  hat.apply.ll <- npplreghat(bws = bw.ll, txdat = tx, tzdat = tz, exdat = ex, ezdat = ez, y = y, output = "apply")
  hat.apply.lp <- npplreghat(bws = bw.lp, txdat = tx, tzdat = tz, exdat = ex, ezdat = ez, y = y, output = "apply")

  expect_equal(as.vector(fit.ll$mean), as.vector(fit.lp$mean), tolerance = 1e-10)
  expect_equal(as.vector(hat.apply.ll), as.vector(hat.apply.lp), tolerance = 1e-10)
  expect_equal(as.vector(H.ll %*% y), as.vector(hat.apply.ll), tolerance = 1e-7)
  expect_equal(as.vector(H.lp %*% y), as.vector(hat.apply.lp), tolerance = 1e-7)
  expect_equal(as.vector(H.ll %*% y), as.vector(fit.ll$mean), tolerance = 1e-7)
  expect_equal(as.vector(H.lp %*% y), as.vector(fit.lp$mean), tolerance = 1e-7)
})
