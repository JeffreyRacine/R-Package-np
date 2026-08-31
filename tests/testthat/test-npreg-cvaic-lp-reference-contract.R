.cvaic_penalty_next_up <- function(value) {
  if (value == 0)
    return(.Machine$double.xmin * .Machine$double.eps)
  bytes <- writeBin(as.double(value), raw(), size = 8L,
                    endian = .Platform$endian)
  if (.Platform$endian != "little")
    bytes <- rev(bytes)
  if (value > 0) {
    carry <- 1L
    for (i in seq_along(bytes)) {
      if (!carry)
        break
      updated <- as.integer(bytes[[i]]) + carry
      bytes[[i]] <- as.raw(updated %% 256L)
      carry <- updated >= 256L
    }
  } else {
    borrow <- 1L
    for (i in seq_along(bytes)) {
      if (!borrow)
        break
      current <- as.integer(bytes[[i]])
      if (current == 0L) {
        bytes[[i]] <- as.raw(255L)
      } else {
        bytes[[i]] <- as.raw(current - 1L)
        borrow <- 0L
      }
    }
  }
  if (.Platform$endian != "little")
    bytes <- rev(bytes)
  readBin(bytes, what = "double", n = 1L, size = 8L,
          endian = .Platform$endian)
}

.cvaic_penalty_c_mse <- function(response) {
  total <- 0
  for (value in response)
    total <- total + value
  center <- total / length(response)
  result <- 0
  for (value in response) {
    difference <- value - center
    result <- result + difference * difference
  }
  result / length(response)
}

.cvaic_penalty_reference <- function(response) {
  n <- length(response)
  if (n <= 3L || any(!is.finite(response)))
    return(NA_real_)
  mse <- .cvaic_penalty_c_mse(response)
  if (!is.finite(mse) || mse <= 0)
    return(NA_real_)
  reference <- log(mse) + (1 + 1 / n) / (1 - 3 / n)
  if (is.finite(reference)) reference else NA_real_
}

.cvaic_penalty_expected <- function(response, multiplier) {
  reference <- .cvaic_penalty_reference(response)
  if (!is.finite(reference))
    return(.Machine$double.xmax)
  candidate <- reference + log(multiplier)
  if (candidate == reference)
    candidate <- .cvaic_penalty_next_up(reference)
  if (!is.finite(candidate) || candidate <= reference ||
      candidate >= .Machine$double.xmax)
    .Machine$double.xmax
  else
    candidate
}

.cvaic_penalty_bad_bandwidth <- function(x, response) {
  npregbw(
    xdat = x, ydat = response,
    bws = 1e-12, bandwidth.compute = FALSE,
    bwtype = "fixed", bwscaling = FALSE,
    ckertype = "epanechnikov", ckerorder = 2L,
    bwmethod = "cv.aic", regtype = "lc"
  )
}

.cvaic_penalty_observe <- function(response, multiplier = 10,
                                    mode = "baseline") {
  x <- data.frame(x = seq(0, 1, length.out = length(response)))
  as.double(np:::.npregbw_eval_only(
    xdat = x, ydat = response,
    bws = .cvaic_penalty_bad_bandwidth(x, response),
    invalid.penalty = mode,
    penalty.multiplier = multiplier
  )$objective)
}

test_that("LP CVAIC agrees with the independent exact-hat definition", {
  old_options <- options(np.tree = FALSE, np.messages = FALSE)
  on.exit(options(old_options), add = TRUE)

  set.seed(202607245L)
  n <- 193L
  x <- data.frame(x1 = runif(n), x2 = runif(n))
  signal <- sin(2 * pi * x$x1) * cos(2 * pi * x$x2) + 0.25 * x$x1 * x$x2
  y <- signal + rnorm(n, sd = 0.35 * sd(signal))
  evaluate_objective <- getFromNamespace(".npregbw_eval_only", "np")

  cases <- list(
    list(type = "fixed", bws = c(0.22, 0.26), degree = c(1L, 1L),
         basis = "glp", bernstein = FALSE),
    list(type = "fixed", bws = c(0.28, 0.32), degree = c(2L, 2L),
         basis = "glp", bernstein = TRUE),
    list(type = "fixed", bws = c(0.30, 0.34), degree = c(2L, 3L),
         basis = "tensor", bernstein = FALSE),
    list(type = "generalized_nn", bws = c(65, 75), degree = c(1L, 1L),
         basis = "glp", bernstein = FALSE),
    list(type = "generalized_nn", bws = c(75, 85), degree = c(2L, 2L),
         basis = "glp", bernstein = FALSE),
    list(type = "adaptive_nn", bws = c(75, 85), degree = c(2L, 2L),
         basis = "glp", bernstein = FALSE)
  )

  for (case in cases) {
    bw <- npregbw(
      xdat = x,
      ydat = y,
      regtype = "lp",
      bwmethod = "cv.aic",
      bwtype = case$type,
      ckertype = "gaussian",
      ckerorder = 2L,
      bws = case$bws,
      bandwidth.compute = FALSE,
      degree = case$degree,
      degree.select = "manual",
      basis = case$basis,
      bernstein.basis = case$bernstein
    )

    native <- evaluate_objective(x, y, bw, objective = "ls")$objective
    hat <- npreghat(bws = bw, txdat = x, output = "matrix")
    fitted <- drop(hat %*% y)
    trace_hat <- sum(diag(hat))
    reference <- log(mean((y - fitted)^2)) +
      (1.0 + trace_hat / n) / (1.0 - (trace_hat + 2.0) / n)

    scale <- max(1.0, abs(reference))
    expect_lt(
      abs(native - reference),
      1e-8 * scale,
      label = paste(
        "type", case$type,
        "degree", paste(case$degree, collapse = ","),
        "basis", case$basis,
        "bernstein", case$bernstein
      )
    )
  }
})

test_that("CVAIC invalid guidance uses the trace-one intercept null", {
  old_options <- options(np.tree = FALSE, np.messages = FALSE)
  on.exit(options(old_options), add = TRUE)

  fixtures <- list(
    n4 = c(-1, -0.25, 0.5, 1.5),
    negative = seq(-0.03, 0.03, length.out = 12L),
    positive = seq(-4, 4, length.out = 12L),
    near_constant = c(rep(1, 11L), 1 + 2^-30)
  )

  for (name in names(fixtures)) {
    response <- fixtures[[name]]
    n <- length(response)
    hat <- matrix(1 / n, nrow = n, ncol = n)
    fitted <- drop(hat %*% response)
    direct <- log(mean((response - fitted)^2)) +
      (1 + sum(diag(hat)) / n) /
      (1 - (sum(diag(hat)) + 2) / n)
    reference <- .cvaic_penalty_reference(response)
    expect_equal(reference, direct, tolerance = 2e-13,
                 info = paste(name, "explicit-hat reference"))

    for (multiplier in c(1, 10)) {
      observed <- .cvaic_penalty_observe(response, multiplier)
      expect_equal(
        observed,
        .cvaic_penalty_expected(response, multiplier),
        tolerance = 8e-15,
        info = paste(name, multiplier)
      )
      expect_gt(observed, reference)
    }

    expect_identical(
      .cvaic_penalty_observe(response, mode = "dbmax"),
      .Machine$double.xmax,
      info = paste(name, "raw-invalid dbmax")
    )
  }
})

test_that("CVAIC terminal, transformation, and state contracts are stable", {
  old_options <- options(np.tree = FALSE, np.messages = FALSE)
  on.exit(options(old_options), add = TRUE)

  x <- data.frame(x = seq(0, 1, length.out = 12L))
  y <- numeric(nrow(x))
  evaluate_objective <- getFromNamespace(".npregbw_eval_only", "np")
  engines <- list(
    list(regtype = "lc", degree = NULL),
    list(regtype = "ll", degree = NULL),
    list(regtype = "lp", degree = 2L)
  )

  for (engine in engines) {
    call <- list(
      xdat = x,
      ydat = y,
      regtype = engine$regtype,
      bwmethod = "cv.aic",
      bwtype = "fixed",
      ckertype = "gaussian",
      ckerorder = 2L,
      bws = 0.3,
      bandwidth.compute = FALSE
    )
    if (!is.null(engine$degree)) call$degree <- engine$degree
    bw <- do.call(npregbw, call)
    baseline <- evaluate_objective(
      x, y, bw, invalid.penalty = "baseline", objective = "ls"
    )$objective
    dbmax <- evaluate_objective(
      x, y, bw, invalid.penalty = "dbmax", objective = "ls"
    )$objective

    expect_identical(as.numeric(baseline), .Machine$double.xmax)
    expect_identical(as.numeric(dbmax), .Machine$double.xmax)
  }

  response <- seq(-0.03, 0.03, length.out = 12L)
  base <- .cvaic_penalty_observe(response)
  translated <- .cvaic_penalty_observe(response + 7)
  scaled <- .cvaic_penalty_observe(-7 * response)
  permuted <- .cvaic_penalty_observe(rev(response))
  expect_equal(translated, base, tolerance = 2e-13)
  expect_equal(scaled, base + 2 * log(7), tolerance = 2e-13)
  expect_equal(permuted, base, tolerance = 2e-14)

  invisible(.cvaic_penalty_observe(response, mode = "dbmax"))
  invisible(.cvaic_penalty_observe(rep(2, length(response))))

  cvls_bw <- npregbw(
    xdat = x, ydat = response,
    bws = 0, bandwidth.compute = FALSE,
    bwtype = "fixed", bwscaling = FALSE,
    ckertype = "uniform", bwmethod = "cv.ls", regtype = "lc"
  )
  invisible(evaluate_objective(
    x, response, cvls_bw,
    invalid.penalty = "baseline", penalty.multiplier = 10
  ))

  binary <- as.double(seq_len(nrow(x)) %% 2L)
  cvks_bw <- npregbw(
    xdat = x, ydat = binary,
    bws = 0, bandwidth.compute = FALSE,
    bwtype = "fixed", bwscaling = FALSE,
    ckertype = "uniform", bwmethod = "cv.ls", regtype = "lc"
  )
  invisible(evaluate_objective(
    x, binary, cvks_bw,
    invalid.penalty = "baseline", penalty.multiplier = 10,
    objective = "ks"
  ))

  expect_identical(.cvaic_penalty_observe(response), base)
})
