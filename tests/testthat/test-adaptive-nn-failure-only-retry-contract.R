adaptive_retry_fixture <- function() {
  n <- 120L
  x1 <- rep(seq(0.04, 0.96, length.out = 6L), each = 20L)
  x2 <- rep(seq(0.06, 0.94, length.out = 5L), times = 24L)
  x3 <- rep(seq(0.08, 0.92, length.out = 4L), each = 5L, times = 6L)
  set.seed(2026090102L)
  list(
    x = data.frame(x1 = x1, x2 = x2, x3 = x3),
    y = sin(2 * pi * x1) + 0.3 * x2 - 0.2 * x3 + rnorm(n, sd = 0.08)
  )
}

test_that("automatic adaptive-NN regression retries only after raw failure", {
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  dat <- adaptive_retry_fixture()

  bw <- npregbw(
    xdat = dat$x, ydat = dat$y,
    bwtype = "adaptive_nn", regtype = "ll",
    nmulti = 1L, itmax = 120L, powell.remin = FALSE
  )
  raw <- getFromNamespace(".npregbw_eval_only", "np")(
    xdat = dat$x, ydat = dat$y, bws = bw,
    invalid.penalty = "dbmax"
  )$objective

  expect_identical(as.double(bw$fval), as.double(raw))
  expect_true(all(is.finite(bw$bw)))
})

test_that("failure-only retry never rewrites explicit adaptive-NN starts", {
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  dat <- adaptive_retry_fixture()

  expect_error(
    npregbw(
      xdat = dat$x, ydat = dat$y, bws = rep.int(12, 3L),
      bwtype = "adaptive_nn", regtype = "ll",
      nmulti = 1L, itmax = 120L, powell.remin = FALSE
    ),
    "optimizer returned a fixed-bandwidth candidate with invalid raw objective",
    fixed = TRUE
  )
})
