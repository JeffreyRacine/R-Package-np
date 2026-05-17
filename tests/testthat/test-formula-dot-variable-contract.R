test_that("formula variables may contain dots without triggering wildcard expansion", {
  set.seed(20260624L)
  old_opts <- options(np.messages = FALSE)
  on.exit(options(old_opts), add = TRUE)

  xdat <- data.frame(x = factor(rbinom(80L, 1L, 0.5)))
  y.irr <- ordered(rbinom(80L, 1L, 0.5))

  bw.dens <- npcdensbw(y.irr ~ x, data = xdat, bws = c(0.2, 0.2),
                       bandwidth.compute = FALSE)
  bw.dist <- npcdistbw(y.irr ~ x, data = xdat, bws = c(0.2, 0.2),
                       bandwidth.compute = FALSE)

  expect_equal(bw.dens$variableNames$response, "y.irr")
  expect_equal(bw.dist$variableNames$response, "y.irr")
  expect_equal(bw.dens$variableNames$terms, "x")
  expect_equal(bw.dist$variableNames$terms, "x")
})

test_that("conditional density formula dot expands only as wildcard", {
  set.seed(20260625L)
  old_opts <- options(np.messages = FALSE)
  on.exit(options(old_opts), add = TRUE)

  dat <- data.frame(
    y = ordered(rbinom(80L, 1L, 0.5)),
    x = factor(rbinom(80L, 1L, 0.5)),
    z = factor(rbinom(80L, 1L, 0.5))
  )

  bw <- npcdensbw(y ~ ., data = dat, bws = c(0.2, 0.2, 0.2),
                  bandwidth.compute = FALSE)

  expect_equal(bw$variableNames$response, "y")
  expect_equal(bw$variableNames$terms, c("x", "z"))
})

test_that("conditional distribution formula dot expands only as wildcard", {
  set.seed(20260626L)
  old_opts <- options(np.messages = FALSE)
  on.exit(options(old_opts), add = TRUE)

  dat <- data.frame(
    y = ordered(rbinom(80L, 1L, 0.5)),
    x = factor(rbinom(80L, 1L, 0.5)),
    z = factor(rbinom(80L, 1L, 0.5))
  )

  bw <- npcdistbw(y ~ ., data = dat, bws = c(0.2, 0.2, 0.2),
                  bandwidth.compute = FALSE)

  expect_equal(bw$variableNames$response, "y")
  expect_equal(bw$variableNames$terms, c("x", "z"))
})
