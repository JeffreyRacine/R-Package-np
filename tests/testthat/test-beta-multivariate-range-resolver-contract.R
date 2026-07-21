.beta_multivariate_training <- function() {
  data.frame(
    x1 = c(0, 0, 0.20, 0.35, 0.60, 0.82, 1, 1),
    x2 = c(-2, -1.5, -1.3, -0.2, 0.5, 1.4, 3, 3)
  )
}

.beta_multivariate_bounds <- function(x) {
  lower <- vapply(x, function(z) {
    u <- sort(unique(as.double(z)))
    u[[1L]] - (u[[2L]] - u[[1L]]) / 2
  }, numeric(1L))
  upper <- vapply(x, function(z) {
    u <- sort(unique(as.double(z)))
    u[[length(u)]] + (u[[length(u)]] - u[[length(u) - 1L]]) / 2
  }, numeric(1L))
  list(lower = lower, upper = upper)
}

test_that("beta range resolves multiple coordinates independently", {
  training <- .beta_multivariate_training()
  expected <- .beta_multivariate_bounds(training)
  dati <- np:::untangle(training)

  resolved <- np:::npKernelBoundsResolve(
    dati = dati,
    varnames = names(training),
    kerbound = "range",
    range.policy = "beta_half_spacing"
  )
  expect_identical(unname(resolved$lb[dati$icon]),
                   unname(expected$lower))
  expect_identical(unname(resolved$ub[dati$icon]),
                   unname(expected$upper))

  range.bw <- npudensbw(
    dat = training, bws = c(0.4, 0.6), bandwidth.compute = FALSE,
    ckertype = "beta", ckerbound = "range"
  )
  fixed.bw <- npudensbw(
    dat = training, bws = c(0.4, 0.6), bandwidth.compute = FALSE,
    ckertype = "beta", ckerbound = "fixed",
    ckerlb = expected$lower, ckerub = expected$upper
  )
  expect_identical(unname(range.bw$ckerlb[range.bw$icon]),
                   unname(expected$lower))
  expect_identical(unname(range.bw$ckerub[range.bw$icon]),
                   unname(expected$upper))

  evaluation <- training[c(1L, 3L, 5L, 7L, 8L), , drop = FALSE]
  range.fit <- npudens(bws = range.bw, tdat = training, edat = evaluation)
  fixed.fit <- npudens(bws = fixed.bw, tdat = training, edat = evaluation)
  expect_identical(fitted(range.fit), fitted(fixed.fit))
  expect_identical(se(range.fit), se(fixed.fit))

  reversed <- training[2:1]
  reversed.expected <- .beta_multivariate_bounds(reversed)
  reversed.bw <- npudensbw(
    dat = reversed, bws = c(0.6, 0.4), bandwidth.compute = FALSE,
    ckertype = "beta", ckerbound = "range"
  )
  expect_identical(unname(reversed.bw$ckerlb[reversed.bw$icon]),
                   unname(reversed.expected$lower))
  expect_identical(unname(reversed.bw$ckerub[reversed.bw$icon]),
                   unname(reversed.expected$upper))
})

test_that("beta range separates data and metadata failures by coordinate", {
  training <- .beta_multivariate_training()

  second.constant <- training
  second.constant$x2 <- 1
  expect_error(
    npudensbw(
      dat = second.constant, bws = c(0.4, 0.6),
      bandwidth.compute = FALSE, ckertype = "beta",
      ckerbound = "range"
    ),
    "Violations: x2",
    fixed = TRUE
  )

  first.constant <- training
  first.constant$x1 <- 1
  expect_error(
    npudensbw(
      dat = first.constant, bws = c(0.4, 0.6),
      bandwidth.compute = FALSE, ckertype = "beta",
      ckerbound = "range"
    ),
    "Violations: x1",
    fixed = TRUE
  )

  malformed <- np:::untangle(training)
  malformed$all.min.next <- malformed$all.min.next[1L]
  expect_error(
    np:::npKernelBoundsResolve(
      dati = malformed,
      varnames = names(training),
      kerbound = "range",
      range.policy = "beta_half_spacing"
    ),
    paste0(
      "Internal beta range metadata inconsistency: expected 2 entries ",
      "in 'all.min.next' and 'all.max.prev'; found 1 and 2"
    ),
    fixed = TRUE
  )
})

test_that("multivariate beta range matches padded fixed bounds unconditionally", {
  training <- .beta_multivariate_training()
  evaluation <- training[c(1L, 3L, 5L, 7L, 8L), , drop = FALSE]
  expected <- .beta_multivariate_bounds(training)
  response <- sin(2 * pi * training$x1) + 0.25 * training$x2

  for (operator in c("normal", "integral")) {
    range.sum <- npksum(
      txdat = training, exdat = evaluation, bws = c(0.4, 0.6),
      ckertype = "beta", ckerbound = "range", operator = operator,
      return.kernel.weights = TRUE
    )
    fixed.sum <- npksum(
      txdat = training, exdat = evaluation, bws = c(0.4, 0.6),
      ckertype = "beta", ckerbound = "fixed",
      ckerlb = expected$lower, ckerub = expected$upper,
      operator = operator, return.kernel.weights = TRUE
    )
    expect_identical(range.sum$ksum, fixed.sum$ksum)
    expect_identical(range.sum$kw, fixed.sum$kw)
  }

  distribution.range.bw <- npudistbw(
    dat = training, bws = c(0.4, 0.6), bandwidth.compute = FALSE,
    ckertype = "beta", ckerbound = "range"
  )
  distribution.fixed.bw <- npudistbw(
    dat = training, bws = c(0.4, 0.6), bandwidth.compute = FALSE,
    ckertype = "beta", ckerbound = "fixed",
    ckerlb = expected$lower, ckerub = expected$upper
  )
  distribution.range <- npudist(
    bws = distribution.range.bw, tdat = training, edat = evaluation
  )
  distribution.fixed <- npudist(
    bws = distribution.fixed.bw, tdat = training, edat = evaluation
  )
  expect_identical(fitted(distribution.range), fitted(distribution.fixed))
  expect_identical(se(distribution.range), se(distribution.fixed))

  regression.range.bw <- npregbw(
    xdat = training, ydat = response, bws = c(0.4, 0.6),
    bandwidth.compute = FALSE, regtype = "lc",
    ckertype = "beta", ckerbound = "range"
  )
  regression.fixed.bw <- npregbw(
    xdat = training, ydat = response, bws = c(0.4, 0.6),
    bandwidth.compute = FALSE, regtype = "lc",
    ckertype = "beta", ckerbound = "fixed",
    ckerlb = expected$lower, ckerub = expected$upper
  )
  regression.range <- npreg(
    bws = regression.range.bw, txdat = training, tydat = response,
    exdat = evaluation
  )
  regression.fixed <- npreg(
    bws = regression.fixed.bw, txdat = training, tydat = response,
    exdat = evaluation
  )
  expect_identical(fitted(regression.range), fitted(regression.fixed))
  expect_identical(se(regression.range), se(regression.fixed))
})

test_that("multivariate beta range resolves conditional X and Y vectors", {
  compare.conditional <- function(x, y, bws) {
    xbounds <- .beta_multivariate_bounds(x)
    ybounds <- .beta_multivariate_bounds(y)
    range.args <- list(
      xdat = x, ydat = y, bws = bws, bandwidth.compute = FALSE,
      cxkertype = "beta", cykertype = "beta",
      cxkerbound = "range", cykerbound = "range"
    )
    fixed.args <- modifyList(range.args, list(
      cxkerbound = "fixed", cxkerlb = xbounds$lower,
      cxkerub = xbounds$upper, cykerbound = "fixed",
      cykerlb = ybounds$lower, cykerub = ybounds$upper
    ))

    density.range.bw <- do.call(npcdensbw, range.args)
    density.fixed.bw <- do.call(npcdensbw, fixed.args)
    distribution.range.bw <- do.call(npcdistbw, range.args)
    distribution.fixed.bw <- do.call(npcdistbw, fixed.args)

    expect_identical(
      unname(density.range.bw$cxkerlb[density.range.bw$ixcon]),
      unname(xbounds$lower)
    )
    expect_identical(
      unname(density.range.bw$cxkerub[density.range.bw$ixcon]),
      unname(xbounds$upper)
    )
    expect_identical(
      unname(density.range.bw$cykerlb[density.range.bw$iycon]),
      unname(ybounds$lower)
    )
    expect_identical(
      unname(density.range.bw$cykerub[density.range.bw$iycon]),
      unname(ybounds$upper)
    )

    rows <- c(1L, 3L, 5L, 8L)
    density.range <- npcdens(
      bws = density.range.bw, txdat = x, tydat = y,
      exdat = x[rows, , drop = FALSE], eydat = y[rows, , drop = FALSE]
    )
    density.fixed <- npcdens(
      bws = density.fixed.bw, txdat = x, tydat = y,
      exdat = x[rows, , drop = FALSE], eydat = y[rows, , drop = FALSE]
    )
    distribution.range <- npcdist(
      bws = distribution.range.bw, txdat = x, tydat = y,
      exdat = x[rows, , drop = FALSE], eydat = y[rows, , drop = FALSE]
    )
    distribution.fixed <- npcdist(
      bws = distribution.fixed.bw, txdat = x, tydat = y,
      exdat = x[rows, , drop = FALSE], eydat = y[rows, , drop = FALSE]
    )

    expect_identical(fitted(density.range), fitted(density.fixed))
    expect_identical(se(density.range), se(density.fixed))
    expect_identical(fitted(distribution.range), fitted(distribution.fixed))
    expect_identical(se(distribution.range), se(distribution.fixed))
  }

  training <- .beta_multivariate_training()
  response <- sin(2 * pi * training$x1) + 0.25 * training$x2

  compare.conditional(
    x = training,
    y = data.frame(y = response),
    bws = c(0.7, 0.4, 0.6)
  )
  compare.conditional(
    x = training["x1"],
    y = data.frame(y1 = response,
                   y2 = training$x2 + 0.3 * training$x1),
    bws = c(0.7, 0.8, 0.4)
  )
})
