transform_bounds_fixture <- function() {
  n <- 36L
  x <- seq(0.03, 0.97, length.out = n)
  u <- factor(rep(c("a", "b", "c"), length.out = n))
  o <- ordered(rep(1:4, length.out = n))
  y <- sin(2 * pi * x) + c(a = 0, b = 0.25, c = -0.15)[u] +
    0.08 * (as.integer(o) - 2.5)
  data.frame(y, x, u, o)
}

test_that("regression transformed search starts on the external bandwidth scale", {
  dat <- transform_bounds_fixture()
  set.seed(314159)
  bw <- npregbw(
    y ~ x + u + o, data = dat,
    scale.factor.init = 0.55, nmulti = 1L,
    transform.bounds = TRUE, powell.remin = FALSE, itmax = 1L
  )

  expect_true(all(is.finite(bw$bw)))
  expect_true(all(bw$bw > 0))
  expect_true(is.finite(bw$fval))

  replay <- np:::.npregbw_eval_only(dat[c("x", "u", "o")], dat$y, bw)
  expect_equal(replay$objective, bw$fval, tolerance = 64 * .Machine$double.eps)
})

test_that("conditional-distribution later starts use transformed coordinates", {
  dat <- transform_bounds_fixture()
  set.seed(271828)
  bw <- npcdistbw(
    y ~ x + u, data = dat,
    scale.factor.init = 0.55, nmulti = 2L,
    transform.bounds = TRUE, powell.remin = FALSE, itmax = 1L
  )

  expect_true(all(is.finite(unlist(bw$sfactor))))
  expect_true(all(unlist(bw$sfactor) > 0))
  expect_true(is.finite(bw$fval))
  expect_length(bw$fval.history, 2L)
  expect_true(all(is.finite(bw$fval.history)))
  expect_equal(bw$fval, min(bw$fval.history), tolerance = 1e-12)

  replay <- np:::.npcdistbw_eval_only(dat[c("x", "u")], dat$y, bws = bw)
  expect_equal(replay$objective, bw$fval, tolerance = 5e-7)
})
