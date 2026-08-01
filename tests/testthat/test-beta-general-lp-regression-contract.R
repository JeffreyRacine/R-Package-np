test_that("manual beta LL and raw degree-one LP share one fitted path", {
  set.seed(2471)
  n <- 61L
  x <- data.frame(x1 = runif(n, 0.03, 0.97), x2 = runif(n, 0.04, 0.96))
  y <- sin(2 * pi * x$x1) + 0.4 * x$x2
  evaluation <- x[c(2L, 11L, 29L, 47L, 59L), , drop = FALSE]
  common <- list(
    xdat = x, ydat = y, bws = c(0.28, 0.31),
    bandwidth.compute = FALSE, bwscaling = FALSE,
    ckertype = "beta", ckerorder = 4L,
    ckerbound = "fixed", ckerlb = c(0, 0), ckerub = c(1, 1)
  )
  ll_bw <- do.call(npregbw, c(common, list(regtype = "ll")))
  lp_bw <- do.call(npregbw, c(
    common,
    list(regtype = "lp", degree = c(1L, 1L),
         bernstein.basis = FALSE)
  ))
  ll <- npreg(
    bws = ll_bw, txdat = x, tydat = y, exdat = evaluation,
    gradients = TRUE
  )
  lp <- npreg(
    bws = lp_bw, txdat = x, tydat = y, exdat = evaluation,
    gradients = TRUE
  )

  expect_identical(fitted(ll), fitted(lp))
  expect_identical(se(ll), se(lp))
  expect_identical(gradients(ll), gradients(lp))
  expect_identical(gradients(ll, errors = TRUE),
                   gradients(lp, errors = TRUE))
})

test_that("beta general LP agrees with independent uncentred WLS", {
  set.seed(3819)
  n <- 83L
  x <- data.frame(x1 = runif(n, 0.04, 0.96), x2 = runif(n, 0.05, 0.95))
  y <- cos(2 * pi * x$x1) + 0.3 * x$x2 + 0.1 * x$x1 * x$x2
  evaluation <- x[c(3L, 17L, 42L, 68L, 79L), , drop = FALSE]
  degree <- c(2L, 1L)
  bw <- npregbw(
    xdat = x, ydat = y, bws = c(0.33, 0.36),
    bandwidth.compute = FALSE, regtype = "lp", basis = "glp",
    degree = degree, bernstein.basis = TRUE, bwscaling = FALSE,
    ckertype = "beta", ckerorder = 6L,
    ckerbound = "fixed", ckerlb = c(0, 0), ckerub = c(1, 1)
  )
  fit <- npreg(bws = bw, txdat = x, tydat = y, exdat = evaluation)
  design <- np:::W.lp(
    xdat = x, degree = degree, basis = "glp", bernstein.basis = TRUE
  )
  evaluation_design <- np:::W.lp(
    xdat = x, exdat = evaluation, degree = degree, basis = "glp",
    bernstein.basis = TRUE
  )
  sums <- npksum(
    bws = bw, txdat = x, exdat = evaluation,
    bandwidth.divide = FALSE, return.kernel.weights = TRUE
  )
  weights <- sums$kw
  if (!is.matrix(weights))
    weights <- matrix(weights, nrow = n, ncol = nrow(evaluation))
  oracle <- vapply(seq_len(nrow(evaluation)), function(index) {
    weight <- weights[, index]
    coefficient <- solve(
      crossprod(design, design * weight),
      crossprod(design, y * weight)
    )
    sum(evaluation_design[index, ] * coefficient)
  }, numeric(1))

  expect_equal(fitted(fit), oracle, tolerance = 2e-9)
})

test_that("automatic beta LP search remains fail-closed", {
  set.seed(5113)
  x <- data.frame(x1 = runif(31), x2 = runif(31))
  y <- x$x1 - x$x2

  expect_error(
    npregbw(
      xdat = x, ydat = y, regtype = "lp", degree = c(1L, 1L),
      ckertype = "beta", ckerorder = 2L,
      ckerbound = "fixed", ckerlb = c(0, 0), ckerub = c(1, 1)
    ),
    "currently supports only regtype = \"lc\""
  )
})
