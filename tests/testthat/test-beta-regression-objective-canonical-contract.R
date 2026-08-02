beta_regression_objective_oracle <- function(training,
                                             response,
                                             bandwidth,
                                             bwtype,
                                             order,
                                             method) {
  weights <- npksum(
    bws = bandwidth,
    txdat = training,
    exdat = training,
    bwtype = bwtype,
    return.kernel.weights = TRUE,
    ckertype = "beta",
    ckerorder = order,
    ckerbound = "fixed",
    ckerlb = rep(0, ncol(training)),
    ckerub = rep(1, ncol(training))
  )$kw

  if (identical(method, "cv.ls")) {
    diag(weights) <- 0
    denominator <- colSums(weights)
    fitted <- colSums(weights * response) / denominator
    return(mean((response - fitted)^2))
  }

  denominator <- colSums(weights)
  fitted <- colSums(weights * response) / denominator
  loss <- mean((response - fitted)^2)
  trace_hat <- sum(diag(weights) / denominator)
  penalty <- (1 + trace_hat / length(response)) /
    (1 - (trace_hat + 2) / length(response))
  if (!is.finite(loss) || loss <= 0 || !is.finite(penalty) || penalty < 0)
    return(.Machine$double.xmax)
  log(loss) + penalty
}

test_that("continuous scalar beta regression objectives use canonical rows", {
  training <- data.frame(
    x1 = c(0.01, 0.035, 0.08, 0.16, 0.29, 0.43,
           0.57, 0.7, 0.82, 0.91, 0.965, 0.995),
    x2 = c(0.02, 0.07, 0.13, 0.21, 0.34, 0.49,
           0.61, 0.73, 0.84, 0.9, 0.96, 0.99)
  )
  response <- sin(2 * pi * training$x1) * cos(pi * training$x2) +
    0.15 * training$x1

  for (order in c(2L, 4L, 6L, 8L)) {
    for (bwtype in c("fixed", "generalized_nn", "adaptive_nn")) {
      bandwidth <- if (identical(bwtype, "fixed")) c(0.17, 0.19) else c(4, 4)
      for (method in c("cv.ls", "cv.aic")) {
        bw <- npregbw(
          xdat = training,
          ydat = response,
          bws = bandwidth,
          bandwidth.compute = FALSE,
          regtype = "lc",
          bwmethod = method,
          bwtype = bwtype,
          bwscaling = FALSE,
          ckertype = "beta",
          ckerorder = order,
          ckerbound = "fixed",
          ckerlb = c(0, 0),
          ckerub = c(1, 1)
        )
        observed <- as.numeric(np:::.npregbw_eval_only(
          xdat = training,
          ydat = response,
          bws = bw,
          invalid.penalty = "dbmax"
        )$objective[[1L]])
        expected <- beta_regression_objective_oracle(
          training, response, bandwidth, bwtype, order, method
        )

        if (identical(expected, .Machine$double.xmax)) {
          expect_identical(observed, expected)
        } else {
          expect_equal(
            observed,
            expected,
            tolerance = 3e-12,
            info = paste(order, bwtype, method)
          )
        }
      }
    }
  }
})

test_that("continuous scalar beta regression activation remains policy-scoped", {
  training <- data.frame(x = seq(0.05, 0.95, length.out = 12L))
  response <- sin(2 * pi * training$x)

  expect_error(
    npregbw(
      xdat = data.frame(x = training$x,
                        group = factor(rep(letters[1:2], each = 6L))),
      ydat = response,
      regtype = "lc",
      ckertype = "beta",
      ckerbound = "fixed",
      ckerlb = 0,
      ckerub = 1
    ),
    "continuous variables only",
    fixed = TRUE
  )
  expect_error(
    npregbw(
      xdat = training,
      ydat = response,
      regtype = "lp",
      degree = 1L,
      ckertype = "beta",
      ckerbound = "fixed",
      ckerlb = 0,
      ckerub = 1
    ),
    "beta regression currently supports only regtype = \"lc\"",
    fixed = TRUE
  )
})
