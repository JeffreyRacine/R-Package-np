library(np)

make_adaptive_nn_conditional_data <- function() {
  set.seed(20260508L)
  n <- 120L
  x <- sort(runif(n))
  y <- sin(5 * x) + 0.3 * x^2 + rnorm(n, sd = 0.08)
  list(
    txdat = data.frame(x = x),
    tydat = data.frame(y = y),
    exdat = data.frame(x = seq(0.1, 0.9, length.out = 40L)),
    eydat = data.frame(y = sin(5 * seq(0.1, 0.9, length.out = 40L)) +
      0.3 * seq(0.1, 0.9, length.out = 40L)^2)
  )
}

make_adaptive_nn_conditional_bw <- function(family, dat, regtype, degree = NULL) {
  fun <- switch(family,
    dens = npcdensbw,
    dist = npcdistbw
  )
  args <- list(
    xdat = dat$txdat,
    ydat = dat$tydat,
    bws = c(12, 12),
    bandwidth.compute = FALSE,
    bwtype = "adaptive_nn",
    regtype = regtype
  )
  if (!is.null(degree)) {
    args$degree <- degree
    args$basis <- "glp"
  }
  do.call(fun, args)
}

eval_adaptive_nn_conditional <- function(family, bw, dat) {
  fun <- switch(family,
    dens = npcdens,
    dist = npcdist
  )
  do.call(fun, list(
    bws = bw,
    txdat = dat$txdat,
    tydat = dat$tydat,
    exdat = dat$exdat,
    eydat = dat$eydat,
    gradients = TRUE
  ))
}

expect_adaptive_nn_lp_contract <- function(family) {
  dat <- make_adaptive_nn_conditional_data()
  fit.ll <- eval_adaptive_nn_conditional(
    family,
    make_adaptive_nn_conditional_bw(family, dat, "ll"),
    dat
  )
  fit.lp1 <- eval_adaptive_nn_conditional(
    family,
    make_adaptive_nn_conditional_bw(family, dat, "lp", 1L),
    dat
  )
  fit.lp2 <- eval_adaptive_nn_conditional(
    family,
    make_adaptive_nn_conditional_bw(family, dat, "lp", 2L),
    dat
  )

  expect_true(all(is.finite(fitted(fit.ll))))
  expect_true(all(is.finite(fitted(fit.lp1))))
  expect_true(all(is.finite(fitted(fit.lp2))))
  expect_equal(fitted(fit.ll), fitted(fit.lp1), tolerance = 1e-10)
  expect_equal(fit.ll$congrad, fit.lp1$congrad, tolerance = 1e-10)

  expect_gt(max(abs(fitted(fit.lp2) - fitted(fit.lp1))), 1e-4)
  expect_gt(max(abs(fit.lp2$congrad - fit.lp1$congrad)), 1e-4)
}

test_that("adaptive-nn conditional distribution lp propagates into npqreg", {
  dat <- make_adaptive_nn_conditional_data()
  q.ll <- npqreg(
    bws = make_adaptive_nn_conditional_bw("dist", dat, "ll"),
    txdat = dat$txdat,
    tydat = dat$tydat,
    exdat = dat$exdat,
    tau = 0.5,
    gradients = TRUE
  )
  q.lp1 <- npqreg(
    bws = make_adaptive_nn_conditional_bw("dist", dat, "lp", 1L),
    txdat = dat$txdat,
    tydat = dat$tydat,
    exdat = dat$exdat,
    tau = 0.5,
    gradients = TRUE
  )
  q.lp2 <- npqreg(
    bws = make_adaptive_nn_conditional_bw("dist", dat, "lp", 2L),
    txdat = dat$txdat,
    tydat = dat$tydat,
    exdat = dat$exdat,
    tau = 0.5,
    gradients = TRUE
  )

  expect_equal(q.ll$quantile, q.lp1$quantile, tolerance = 1e-10)
  expect_equal(q.ll$quanterr, q.lp1$quanterr, tolerance = 1e-10)
  expect_equal(q.ll$quantgrad, q.lp1$quantgrad, tolerance = 1e-10)
  expect_equal(q.ll$quantgerr, q.lp1$quantgerr, tolerance = 1e-10)

  expect_gt(max(abs(q.lp2$quantile - q.lp1$quantile)), 1e-4)
  expect_gt(max(abs(q.lp2$quanterr - q.lp1$quanterr)), 1e-4)
  expect_gt(max(abs(q.lp2$quantgrad - q.lp1$quantgrad)), 1e-4)
  expect_gt(max(abs(q.lp2$quantgerr - q.lp1$quantgerr)), 1e-4)
})

test_that("adaptive-nn conditional density lp honors degree metadata", {
  expect_adaptive_nn_lp_contract("dens")
})

test_that("adaptive-nn conditional distribution lp honors degree metadata", {
  expect_adaptive_nn_lp_contract("dist")
})
