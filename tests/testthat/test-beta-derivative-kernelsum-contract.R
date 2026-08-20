beta_pdf_derivative_oracle <- function(x, u, h, order, a = 0, b = 1) {
  L <- b - a
  xt <- (x - a) / L
  ut <- (u - a) / L
  coefficients <- switch(as.character(order),
                         `2` = 1,
                         `4` = c(2, -1),
                         `6` = c(3, -3, 1),
                         `8` = c(4, -6, 4, -1))
  tau <- (L / h)^2
  result <- numeric(length(u))
  for (scale in seq_along(coefficients)) {
    concentration <- tau / scale
    alpha <- 1 + xt * concentration
    beta <- 1 + (1 - xt) * concentration
    value <- dbeta(ut, alpha, beta) / L
    score <- concentration / L *
      (log(ut) - log1p(-ut) - digamma(alpha) + digamma(beta))
    result <- result + coefficients[scale] * value * score
  }
  result
}

test_that("beta PDF derivatives match analytic oracles at every order", {
  training <- data.frame(x = c(.017, .08, .21, .47, .73, .94, .989))
  evaluation <- data.frame(x = c(.003, .11, .37, .5, .82, .997))

  for (order in c(2L, 4L, 6L, 8L)) {
    fit <- npksum(
      bws = .14, txdat = training, exdat = evaluation,
      operator = "derivative", return.kernel.weights = TRUE,
      ckertype = "beta", ckerorder = order,
      ckerbound = "fixed", ckerlb = 0, ckerub = 1
    )
    expected <- vapply(
      evaluation$x, beta_pdf_derivative_oracle,
      numeric(nrow(training)), u = training$x, h = .14, order = order
    )
    expect_equal(fit$kw, expected, tolerance = 2e-11)
    expect_equal(as.double(fit$ksum), colSums(expected), tolerance = 2e-11)
  }
})

test_that("beta direct and permutation derivatives agree", {
  training <- data.frame(
    x = c(.017, .08, .21, .47, .73, .94, .989),
    z = c(.9, .18, .67, .31, .79, .42, .06)
  )
  evaluation <- data.frame(x = c(.11, .37, .82), z = c(.2, .7, .44))
  common <- list(
    bws = c(.14, .18), txdat = training, exdat = evaluation,
    ckertype = "beta", ckerorder = 6,
    ckerbound = "fixed", ckerlb = c(0, 0), ckerub = c(1, 1)
  )
  permutation <- do.call(npksum, c(common, list(
    permutation.operator = "derivative", return.kernel.weights = TRUE,
    return.derivative.kernel.weights = TRUE
  )))
  for (dimension in 1:2) {
    operators <- rep("normal", 2)
    operators[dimension] <- "derivative"
    direct <- do.call(npksum, c(common, list(
      operator = operators, return.kernel.weights = TRUE
    )))
    expect_equal(as.double(permutation$p.ksum[, dimension]),
                 as.double(direct$ksum), tolerance = 2e-12)
    expect_equal(permutation$p.kw[, , dimension], direct$kw,
                 tolerance = 2e-12)
  }
})

test_that("powered beta derivatives match raw-weight oracles", {
  training <- data.frame(
    x = c(.07, .16, .29, .43, .61, .76, .91),
    z = c(.88, .19, .68, .34, .79, .47, .11)
  )
  evaluation <- data.frame(x = c(.12, .31, .52, .81),
                           z = c(.22, .64, .39, .73))
  arguments <- list(
    bws = c(.18, .2), txdat = training, exdat = evaluation,
    operator = c("derivative", "normal"),
    return.kernel.weights = TRUE,
    ckertype = "beta", ckerorder = 6,
    ckerbound = "fixed", ckerlb = c(0, 0), ckerub = c(1, 1)
  )
  raw <- do.call(npksum, arguments)

  for (power in c(-2L, -1L, 0L, 2L, 3L)) {
    powered <- do.call(npksum, c(arguments, list(kernel.pow = power)))
    expected_weights <- ifelse(raw$kw == 0, 0, raw$kw^power)
    expect_equal(as.double(powered$ksum), colSums(expected_weights),
                 tolerance = 3e-10)
    expect_identical(powered$kw, raw$kw)
  }
})

test_that("powered beta derivative tensors retain public layout", {
  training <- data.frame(
    x = c(.07, .16, .29, .43, .61, .76, .91),
    z = c(.88, .19, .68, .34, .79, .47, .11)
  )
  evaluation <- data.frame(x = c(.12, .31, .52, .81),
                           z = c(.22, .64, .39, .73))
  response <- cbind(y1 = seq_len(nrow(training)),
                    y2 = cos(seq_len(nrow(training))))
  case_weights <- cbind(w1 = 1, w2 = seq(.6, 1.2, length.out = nrow(training)))
  arguments <- list(
    bws = c(.18, .2), txdat = training, tydat = response,
    weights = case_weights, exdat = evaluation,
    operator = c("derivative", "normal"),
    return.kernel.weights = TRUE,
    ckertype = "beta", ckerorder = 8,
    ckerbound = "fixed", ckerlb = c(0, 0), ckerub = c(1, 1)
  )
  raw <- do.call(npksum, arguments)
  powered <- do.call(npksum, c(arguments, list(kernel.pow = 2L)))
  powered_weights <- raw$kw^2

  for (evaluation_index in seq_len(nrow(evaluation)))
    for (response_column in seq_len(ncol(response)))
      for (weight_column in seq_len(ncol(case_weights)))
        expect_equal(
          powered$ksum[weight_column, response_column, evaluation_index],
          sum(powered_weights[, evaluation_index] *
              response[, response_column] *
              case_weights[, weight_column]),
          tolerance = 5e-10
        )
  expect_identical(powered$kw, raw$kw)
})

test_that("powered beta permutation derivatives reuse direct rows", {
  training <- data.frame(
    x = c(.07, .16, .29, .43, .61, .76, .91),
    z = c(.88, .19, .68, .34, .79, .47, .11)
  )
  evaluation <- data.frame(x = c(.12, .31, .52, .81),
                           z = c(.22, .64, .39, .73))
  common <- list(
    bws = c(.18, .2), txdat = training, exdat = evaluation,
    kernel.pow = 2L, ckertype = "beta", ckerorder = 8,
    ckerbound = "fixed", ckerlb = c(0, 0), ckerub = c(1, 1)
  )
  permutation <- do.call(npksum, c(common, list(
    permutation.operator = "derivative", return.kernel.weights = TRUE,
    return.derivative.kernel.weights = TRUE
  )))

  for (dimension in 1:2) {
    operators <- rep("normal", 2)
    operators[dimension] <- "derivative"
    direct <- do.call(npksum, c(common, list(
      operator = operators, return.kernel.weights = TRUE
    )))
    expect_equal(as.double(permutation$p.ksum[, dimension]),
                 as.double(direct$ksum), tolerance = 3e-10)
    expect_identical(permutation$p.kw[, , dimension], direct$kw)
  }
})

test_that("powered beta endpoint jumps retain structural semantics", {
  arguments <- list(
    bws = .14, txdat = data.frame(x = c(0, .2, .8, 1)),
    exdat = data.frame(x = c(0, 1)), operator = "derivative",
    return.kernel.weights = TRUE, ckertype = "beta", ckerorder = 8,
    ckerbound = "fixed", ckerlb = 0, ckerub = 1
  )
  raw <- NULL
  expect_warning(raw <- do.call(npksum, arguments), "infinite endpoint")

  squared <- NULL
  expect_warning(squared <- do.call(npksum, c(arguments, list(
    kernel.pow = 2L
  ))), "infinite endpoint")
  expect_true(all(is.infinite(squared$ksum)))
  expect_true(all(squared$ksum > 0))
  expect_identical(squared$kw, raw$kw)

  cubed <- NULL
  expect_warning(cubed <- do.call(npksum, c(arguments, list(
    kernel.pow = 3L
  ))), "infinite endpoint")
  expect_identical(sign(as.double(cubed$ksum)), c(-1, 1))
  expect_identical(cubed$kw, raw$kw)

  zeroth <- do.call(npksum, c(arguments, list(kernel.pow = 0L)))
  expected_zero <- ifelse(raw$kw == 0, 0, 1)
  expect_equal(as.double(zeroth$ksum), colSums(expected_zero),
               tolerance = 2e-12)
  expect_identical(zeroth$kw, raw$kw)
})

beta_derivative_mpi_domains <- function() {
  list(
    mode = c("fixed", "generalized_nn", "adaptive_nn"),
    order = c(2L, 4L, 6L, 8L),
    derivative_dimension = 1:2,
    companion = c("normal", "integral", "convolution"),
    power = c(0L, 2L, 3L)
  )
}

beta_derivative_mpi_pairwise_cases <- function() {
  data.frame(
    mode = c(
      "fixed", "generalized_nn", "adaptive_nn",
      "fixed", "generalized_nn", "adaptive_nn",
      "fixed", "generalized_nn", "adaptive_nn",
      "fixed", "generalized_nn", "adaptive_nn"
    ),
    order = rep(c(2L, 4L, 6L, 8L), each = 3L),
    derivative_dimension = c(2L, 2L, 1L, 2L, 2L, 1L,
                             1L, 2L, 1L, 2L, 1L, 2L),
    companion = c(
      "integral", "normal", "convolution",
      "convolution", "integral", "normal",
      "normal", "convolution", "integral",
      "integral", "normal", "convolution"
    ),
    power = c(3L, 2L, 0L, 2L, 0L, 3L, 0L, 3L, 2L, 3L, 2L, 0L),
    stringsAsFactors = FALSE
  )
}

test_that("MPI powered beta derivative cases cover every route-axis pair", {
  domains <- beta_derivative_mpi_domains()
  cases <- beta_derivative_mpi_pairwise_cases()
  expect_identical(anyDuplicated(cases), 0L)
  for (pair in combn(names(cases), 2L, simplify = FALSE)) {
    observed <- unique(cases[pair])
    complete <- expand.grid(domains[pair], stringsAsFactors = FALSE)
    expect_setequal(do.call(paste, observed), do.call(paste, complete))
  }
})

test_that("powered beta derivatives span orders, operators, and bandwidth modes", {
  training <- data.frame(
    x = c(.07, .16, .29, .43, .61, .76, .91),
    z = c(.88, .19, .68, .34, .79, .47, .11)
  )
  evaluation <- data.frame(x = c(.12, .31, .52, .81),
                           z = c(.22, .64, .39, .73))
  cases <- beta_derivative_mpi_pairwise_cases()
  for (case_index in seq_len(nrow(cases))) {
    mode <- cases$mode[[case_index]]
    order <- cases$order[[case_index]]
    derivative_dimension <- cases$derivative_dimension[[case_index]]
    companion <- cases$companion[[case_index]]
    power <- cases$power[[case_index]]
    bandwidth <- if (identical(mode, "fixed")) c(.18, .2) else c(4, 4)
    operators <- rep(companion, 2)
    operators[derivative_dimension] <- "derivative"
    arguments <- list(
      bws = bandwidth, txdat = training, exdat = evaluation,
      bwtype = mode, operator = operators,
      return.kernel.weights = TRUE, bandwidth.divide = TRUE,
      ckertype = "beta", ckerorder = order,
      ckerbound = "fixed", ckerlb = c(0, 0), ckerub = c(1, 1)
    )
    raw <- do.call(npksum, arguments)
    powered <- do.call(npksum, c(arguments, list(kernel.pow = power)))
    expected_weights <- ifelse(raw$kw == 0, 0, raw$kw^power)
    expect_equal(as.double(powered$ksum), colSums(expected_weights),
                 tolerance = 8e-9)
    expect_identical(powered$kw, raw$kw)
  }
})

test_that("powered beta derivative deletion remains structural", {
  training <- data.frame(
    x = c(.07, .16, .29, .43, .61, .76, .91),
    z = c(.88, .19, .68, .34, .79, .47, .11)
  )

  for (mode in c("fixed", "generalized_nn", "adaptive_nn")) {
    bandwidth <- if (identical(mode, "fixed")) c(.18, .2) else c(4, 4)
    arguments <- list(
      bws = bandwidth, txdat = training, bwtype = mode,
      operator = c("derivative", "convolution"),
      leave.one.out = TRUE, return.kernel.weights = TRUE,
      ckertype = "beta", ckerorder = 8,
      ckerbound = "fixed", ckerlb = c(0, 0), ckerub = c(1, 1)
    )
    raw <- do.call(npksum, arguments)
    powered <- do.call(npksum, c(arguments, list(kernel.pow = 2L)))
    expect_true(all(diag(raw$kw) == 0))
    expect_identical(powered$kw, raw$kw)
    expect_equal(as.double(powered$ksum), colSums(raw$kw^2),
                 tolerance = 8e-9)
  }
})

test_that("forced tree mode preserves canonical beta permutation derivatives", {
  training <- data.frame(
    x = c(.017, .08, .21, .47, .73, .94, .989),
    z = c(.9, .18, .67, .31, .79, .42, .06)
  )
  evaluation <- data.frame(x = c(.11, .37, .82), z = c(.2, .7, .44))
  arguments <- list(
    bws = c(.14, .18), txdat = training, exdat = evaluation,
    permutation.operator = "derivative",
    return.kernel.weights = TRUE,
    return.derivative.kernel.weights = TRUE,
    ckertype = "beta", ckerorder = 8,
    ckerbound = "fixed", ckerlb = c(0, 0), ckerub = c(1, 1)
  )
  old_options <- options(np.tree = FALSE)
  on.exit(options(old_options), add = TRUE)
  without_tree <- do.call(npksum, arguments)
  options(np.tree = TRUE)
  with_tree <- do.call(npksum, arguments)

  expect_identical(with_tree$ksum, without_tree$ksum)
  expect_identical(with_tree$kw, without_tree$kw)
  expect_identical(with_tree$p.ksum, without_tree$p.ksum)
  expect_identical(with_tree$p.kw, without_tree$p.kw)
})

test_that("beta scalar endpoint derivative signs are explicit", {
  fit <- NULL
  expect_warning(fit <- npksum(
    bws = .14, txdat = data.frame(x = c(0, .2, .8, 1)),
    exdat = data.frame(x = c(0, 1)), operator = "derivative",
    return.kernel.weights = TRUE, ckertype = "beta", ckerorder = 8,
    ckerbound = "fixed", ckerlb = 0, ckerub = 1
  ), "infinite endpoint")
  expect_identical(sign(fit$kw[cbind(c(1L, 4L), c(1L, 2L))]), c(-1, 1))
  expect_identical(fit$kw[cbind(c(4L, 1L), c(1L, 2L))], c(0, 0))
})

test_that("beta endpoint derivatives without matching observations are one-sided", {
  training <- data.frame(x = c(.07, .19, .38, .61, .83, .94))
  step <- 2e-6

  for (order in c(2L, 4L, 6L, 8L)) {
    evaluate <- function(x, operator = "normal") {
      as.double(npksum(
        bws = .16, txdat = training, exdat = data.frame(x = x),
        operator = operator, ckertype = "beta", ckerorder = order,
        ckerbound = "fixed", ckerlb = 0, ckerub = 1
      )$ksum)
    }
    derivative <- evaluate(c(0, 1), "derivative")
    oracle <- c(
      (evaluate(step) - evaluate(0)) / step,
      (evaluate(1) - evaluate(1 - step)) / step
    )
    expect_true(all(is.finite(derivative)))
    expect_equal(derivative, oracle, tolerance = 3e-4)
  }
})
