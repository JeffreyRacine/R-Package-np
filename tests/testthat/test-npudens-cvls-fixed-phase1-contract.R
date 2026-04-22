library(npRmpi)

bounded_gaussian_kernel <- function(eval_x, train_x, h, lb, ub, order) {
  z <- (train_x - eval_x) / h
  pdf <- switch(
    as.character(as.integer(order)),
    `2` = dnorm(z),
    `4` = (1.5 - 0.5 * z^2) * dnorm(z),
    stop("unsupported order")
  )
  cdf <- switch(
    as.character(as.integer(order)),
    `2` = pnorm,
    `4` = function(v) pnorm(v) + 0.5 * v * dnorm(v),
    stop("unsupported order")
  )
  denom <- h * (cdf((ub - eval_x) / h) - cdf((lb - eval_x) / h))
  pdf / denom
}

reconstruct_npudens_cvls_fixed <- function(x, bw, q = 1201L) {
  stopifnot(
    identical(bw$method, "cv.ls"),
    identical(bw$type, "fixed"),
    identical(bw$ckertype, "gaussian"),
    bw$ncon == 1L,
    bw$nuno == 0L,
    bw$nord == 0L
  )

  h <- as.numeric(bw$bw[1])
  lb <- as.numeric(bw$ckerlb[bw$icon][1])
  ub <- as.numeric(bw$ckerub[bw$icon][1])
  grid <- seq(lb, ub, length.out = q)
  trap <- diff(grid)[1] * c(0.5, rep(1, q - 2L), 0.5)

  fit <- npudens(
    bws = bw,
    tdat = data.frame(x = x),
    edat = data.frame(x = grid)
  )

  i1 <- sum(trap * as.numeric(fit$dens)^2)
  i2 <- mean(vapply(seq_along(x), function(i) {
    mean(bounded_gaussian_kernel(x[i], x[-i], h, lb, ub, bw$ckerorder))
  }, numeric(1)))

  list(i1 = i1, i2 = i2, objective = i1 - 2 * i2)
}

test_that("phase1 npudensbw cv.ls fixed bounded objective matches direct reconstruction", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  set.seed(20260421)
  x <- rbeta(64L, 2, 5)
  dat <- data.frame(x = x)

  bw2 <- npudensbw(
    dat = dat,
    bwmethod = "cv.ls",
    bwtype = "fixed",
    ckertype = "gaussian",
    ckerorder = 2L,
    ckerbound = "range",
    nmulti = 1
  )
  rec2 <- reconstruct_npudens_cvls_fixed(x, bw2)
  expect_true(is.finite(as.numeric(bw2$bw[1])))
  expect_true(as.numeric(bw2$bw[1]) > 0)
  expect_equal(as.numeric(bw2$fval), -rec2$objective, tolerance = 1e-2)

  bw4 <- npudensbw(
    dat = dat,
    bwmethod = "cv.ls",
    bwtype = "fixed",
    ckertype = "gaussian",
    ckerorder = 4L,
    ckerbound = "range",
    nmulti = 1
  )
  rec4 <- reconstruct_npudens_cvls_fixed(x, bw4)
  expect_true(is.finite(as.numeric(bw4$bw[1])))
  expect_true(as.numeric(bw4$bw[1]) > 0)
  expect_equal(as.numeric(bw4$fval), -rec4$objective, tolerance = 1e-2)
})
