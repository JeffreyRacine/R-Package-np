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
