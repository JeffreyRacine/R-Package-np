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

beta_regression_lp_objective_oracle <- function(training,
                                                response,
                                                bandwidth,
                                                bwtype,
                                                order,
                                                method,
                                                degree,
                                                basis,
                                                bernstein) {
  design <- np:::W.lp(
    xdat = training,
    degree = degree,
    basis = basis,
    bernstein.basis = bernstein
  )
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
  n <- nrow(training)
  loss <- 0
  trace_hat <- 0

  for (evaluation in seq_len(n)) {
    weight <- weights[, evaluation]
    gram <- crossprod(design, design * weight)
    projection <- solve(gram, design[evaluation, ])
    influence <- weight * drop(design %*% projection)
    fitted <- sum(influence * response)
    residual <- response[evaluation] - fitted

    if (identical(method, "cv.ls")) {
      residual <- residual / (1 - influence[evaluation])
      loss <- loss + residual^2
    } else {
      loss <- loss + residual^2
      trace_hat <- trace_hat + influence[evaluation]
    }
  }

  loss <- loss / n
  if (identical(method, "cv.ls"))
    return(loss)
  penalty <- (1 + trace_hat / n) / (1 - (trace_hat + 2) / n)
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

test_that("continuous general LP beta objectives use canonical full rows", {
  set.seed(781)
  n <- 67L
  training <- data.frame(
    x1 = runif(n, 0.02, 0.98),
    x2 = runif(n, 0.02, 0.98)
  )
  response <- sin(2 * pi * training$x1) + 0.4 * training$x2 +
    rnorm(n, sd = 0.08)
  bases <- c("glp", "additive", "tensor")

  for (order in c(2L, 4L, 6L, 8L)) {
    for (bwtype in c("fixed", "generalized_nn", "adaptive_nn")) {
      bandwidth <- if (identical(bwtype, "fixed")) c(0.42, 0.44) else c(49, 51)
      for (basis in bases) {
        degree <- if (identical(basis, "tensor")) c(1L, 1L) else c(2L, 1L)
        for (bernstein in c(FALSE, TRUE)) {
          for (method in c("cv.ls", "cv.aic")) {
            bw <- npregbw(
              xdat = training,
              ydat = response,
              bws = bandwidth,
              bandwidth.compute = FALSE,
              regtype = "lp",
              degree = degree,
              basis = basis,
              bernstein.basis = bernstein,
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
            expected <- beta_regression_lp_objective_oracle(
              training, response, bandwidth, bwtype, order, method,
              degree, basis, bernstein
            )
            info <- paste(order, bwtype, basis, bernstein, method)

            if (identical(expected, .Machine$double.xmax)) {
              expect_identical(observed, expected, info = info)
            } else {
              expect_equal(observed, expected, tolerance = 1e-8, info = info)
            }
          }
        }
      }
    }
  }
})

test_that("beta ll and raw degree-one LP objectives share one canonical route", {
  set.seed(811)
  n <- 43L
  training <- data.frame(x1 = runif(n), x2 = runif(n))
  response <- cos(2 * pi * training$x1) + training$x2

  for (bwtype in c("fixed", "generalized_nn", "adaptive_nn")) {
    bandwidth <- if (identical(bwtype, "fixed")) c(0.31, 0.34) else c(29, 31)
    for (method in c("cv.ls", "cv.aic")) {
      common <- list(
        xdat = training,
        ydat = response,
        bws = bandwidth,
        bandwidth.compute = FALSE,
        bwmethod = method,
        bwtype = bwtype,
        bwscaling = FALSE,
        ckertype = "beta",
        ckerorder = 6L,
        ckerbound = "fixed",
        ckerlb = c(0, 0),
        ckerub = c(1, 1)
      )
      ll <- do.call(npregbw, c(common, list(regtype = "ll")))
      lp <- do.call(npregbw, c(common, list(
        regtype = "lp",
        degree = c(1L, 1L),
        basis = "glp",
        bernstein.basis = FALSE
      )))
      ll_objective <- np:::.npregbw_eval_only(
        training, response, ll, invalid.penalty = "dbmax"
      )$objective[[1L]]
      lp_objective <- np:::.npregbw_eval_only(
        training, response, lp, invalid.penalty = "dbmax"
      )$objective[[1L]]

      expect_identical(ll$regtype.engine, "lp")
      expect_identical(ll$basis.engine, "glp")
      expect_identical(ll$degree.engine, c(1L, 1L))
      expect_identical(ll_objective, lp_objective)
    }
  }
})

test_that("continuous beta regression activation remains categorically scoped", {
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
})
