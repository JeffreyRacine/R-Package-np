test_that("beta regression gradients match finite differences and hats", {
  training <- data.frame(
    x = c(.03, .08, .16, .27, .41, .56, .69, .81, .91, .98),
    z = c(.07, .91, .22, .75, .38, .62, .14, .86, .49, .31)
  )
  response <- sin(2.4 * training$x) + .6 * training$z^2
  evaluation <- data.frame(x = c(.13, .34, .62, .87),
                           z = c(.2, .47, .78, .39))
  step <- 2e-6

  for (order in c(2L, 4L, 6L, 8L)) {
    common <- list(
      bws = c(.17, .21), regtype = "lc", ckertype = "beta",
      ckerorder = order, ckerbound = "fixed",
      ckerlb = c(0, 0), ckerub = c(1, 1)
    )
    fit <- do.call(npreg, c(list(
      txdat = training, tydat = response, exdat = evaluation,
      gradients = TRUE
    ), common))
    plus <- minus <- evaluation
    plus$x <- plus$x + step
    minus$x <- minus$x - step
    oracle <- (
      fitted(do.call(npreg, c(list(txdat = training, tydat = response,
                                   exdat = plus), common))) -
      fitted(do.call(npreg, c(list(txdat = training, tydat = response,
                                   exdat = minus), common)))
    ) / (2 * step)
    expect_equal(gradients(fit)[, 1L], oracle, tolerance = 3e-6)

    bw <- do.call(npregbw, c(list(
      xdat = training, ydat = response, bandwidth.compute = FALSE
    ), common))
    H <- npreghat(bws = bw, txdat = training, exdat = evaluation,
                  output = "matrix", s = c(1L, 0L))
    expect_equal(as.double(H %*% response), gradients(fit)[, 1L],
                 tolerance = 2e-9)
  }
})

test_that("beta ratio derivatives honor endpoint structural cancellation", {
  training <- data.frame(x = c(0, .12, .3, .55, .82, 1))
  evaluation <- data.frame(x = c(0, .4, 1))

  constant <- npreg(
    bws = .16, txdat = training, tydat = rep(7, nrow(training)),
    exdat = evaluation, gradients = TRUE, regtype = "lc",
    ckertype = "beta", ckerorder = 8,
    ckerbound = "fixed", ckerlb = 0, ckerub = 1
  )
  expect_identical(as.double(gradients(constant)), c(0, 0, 0))
  expect_identical(as.double(constant$gerr), c(0, 0, 0))

  varying <- NULL
  expect_warning(varying <- npreg(
    bws = .16, txdat = training, tydat = c(4, 1, 2, 3, 5, -1),
    exdat = evaluation, gradients = TRUE, regtype = "lc",
    ckertype = "beta", ckerorder = 8,
    ckerbound = "fixed", ckerlb = 0, ckerub = 1
  ), "infinite endpoint")
  expect_true(is.infinite(gradients(varying)[1L, 1L]))
  expect_true(is.infinite(gradients(varying)[3L, 1L]))
  expect_true(is.finite(gradients(varying)[2L, 1L]))
})

test_that("beta regression endpoint gradients have finite one-sided limits", {
  training <- data.frame(x = c(.06, .17, .31, .48, .66, .84, .95))
  response <- cos(2.2 * training$x) + training$x^2
  common <- list(
    bws = .16, txdat = training, tydat = response, regtype = "lc",
    ckertype = "beta", ckerorder = 8,
    ckerbound = "fixed", ckerlb = 0, ckerub = 1
  )
  step <- 1e-7
  fit <- do.call(npreg, c(common, list(
    exdat = data.frame(x = c(0, 1)), gradients = TRUE
  )))
  evaluate <- function(x) fitted(do.call(npreg, c(common, list(
    exdat = data.frame(x = x)
  ))))
  oracle <- c(
    (evaluate(step) - evaluate(0)) / step,
    (evaluate(1) - evaluate(1 - step)) / step
  )
  actual <- as.double(gradients(fit))
  expect_true(all(is.finite(actual)))
  expect_lt(max(abs(actual - oracle) / pmax(1, abs(actual))), 6e-5)
})

test_that("conditional beta gradients match finite differences", {
  training_x <- data.frame(x = c(.03, .08, .16, .27, .41, .56, .69, .81, .91, .98))
  training_y <- data.frame(y = c(.06, .19, .11, .36, .52, .43, .77, .64, .88, .95))
  evaluation_x <- data.frame(x = c(.13, .34, .62, .87))
  evaluation_y <- data.frame(y = c(.17, .45, .71, .9))
  step <- 2e-6

  for (distribution in c(FALSE, TRUE)) {
    bwfun <- if (distribution) npcdistbw else npcdensbw
    fitfun <- if (distribution) npcdist else npcdens
    bw <- bwfun(
      xdat = training_x, ydat = training_y, bws = c(.16, .18),
      bandwidth.compute = FALSE,
      cxkertype = "beta", cxkerorder = 6,
      cxkerbound = "fixed", cxkerlb = 0, cxkerub = 1,
      cykertype = "beta", cykerorder = 8,
      cykerbound = "fixed", cykerlb = 0, cykerub = 1
    )
    fit <- fitfun(bws = bw, txdat = training_x, tydat = training_y,
                  exdat = evaluation_x, eydat = evaluation_y,
                  gradients = TRUE)
    plus <- minus <- evaluation_x
    plus$x <- plus$x + step
    minus$x <- minus$x - step
    oracle <- (
      fitted(fitfun(bws = bw, txdat = training_x, tydat = training_y,
                    exdat = plus, eydat = evaluation_y)) -
      fitted(fitfun(bws = bw, txdat = training_x, tydat = training_y,
                    exdat = minus, eydat = evaluation_y))
    ) / (2 * step)
    expect_equal(as.double(gradients(fit)), oracle, tolerance = 5e-6)
    expect_true(all(is.finite(fit$congerr) & fit$congerr >= 0))
  }
})

test_that("conditional beta endpoint gradients have finite one-sided limits", {
  training_x <- data.frame(x = c(.06, .14, .27, .39, .53, .68, .81, .94))
  training_y <- data.frame(y = c(.11, .29, .18, .47, .39, .72, .61, .88))
  evaluation_y <- data.frame(y = c(.36, .74))
  step <- 1e-7

  for (distribution in c(FALSE, TRUE)) {
    bwfun <- if (distribution) npcdistbw else npcdensbw
    fitfun <- if (distribution) npcdist else npcdens
    bw <- bwfun(
      xdat = training_x, ydat = training_y, bws = c(.16, .18),
      bandwidth.compute = FALSE,
      cxkertype = "beta", cxkerorder = 8,
      cxkerbound = "fixed", cxkerlb = 0, cxkerub = 1,
      cykertype = "beta", cykerorder = 8,
      cykerbound = "fixed", cykerlb = 0, cykerub = 1
    )
    evaluate <- function(x, y) fitted(fitfun(
      bws = bw, txdat = training_x, tydat = training_y,
      exdat = data.frame(x = x), eydat = data.frame(y = y)
    ))
    fit <- fitfun(
      bws = bw, txdat = training_x, tydat = training_y,
      exdat = data.frame(x = c(0, 1)), eydat = evaluation_y,
      gradients = TRUE
    )
    oracle <- c(
      (evaluate(step, evaluation_y$y[1L]) -
         evaluate(0, evaluation_y$y[1L])) / step,
      (evaluate(1, evaluation_y$y[2L]) -
         evaluate(1 - step, evaluation_y$y[2L])) / step
    )
    actual <- as.double(gradients(fit))
    expect_true(all(is.finite(actual)))
    expect_lt(max(abs(actual - oracle) / pmax(1, abs(actual))), 1e-4)
  }
})

test_that("beta derivative plot-data routes retain native kernel descriptors", {
  x <- data.frame(x = c(.04, .1, .18, .29, .43, .58, .7, .82, .91, .97))
  y <- sin(2.7 * x$x) + .2 * x$x^2
  cy <- data.frame(y = c(.07, .2, .13, .38, .49, .45, .73, .66, .86, .94))
  rbw <- npregbw(
    xdat = x, ydat = y, bws = .16, bandwidth.compute = FALSE,
    regtype = "lc", ckertype = "beta", ckerorder = 8,
    ckerbound = "fixed", ckerlb = 0, ckerub = 1
  )
  regression <- plot(
    rbw, xdat = x, ydat = y, gradients = TRUE, errors = "none",
    output = "data", neval = 11L
  )
  expect_true(is.list(regression) && length(regression) > 0L)

  for (bwfun in list(npcdensbw, npcdistbw)) {
    bw <- bwfun(
      xdat = x, ydat = cy, bws = c(.18, .16),
      bandwidth.compute = FALSE,
      cxkertype = "beta", cxkerorder = 8,
      cxkerbound = "fixed", cxkerlb = 0, cxkerub = 1,
      cykertype = "beta", cykerorder = 8,
      cykerbound = "fixed", cykerlb = 0, cykerub = 1
    )
    conditional <- plot(
      bw, xdat = x, ydat = cy, gradients = TRUE, errors = "none",
      output = "data", neval = 11L
    )
    expect_true(is.list(conditional) && length(conditional) > 0L)
  }
})
