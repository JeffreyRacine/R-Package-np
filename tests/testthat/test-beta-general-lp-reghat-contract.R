beta_lp_hat_oracle <- function(bws, training, evaluation, degree,
                               basis = "glp", bernstein = FALSE,
                               derivative = NULL) {
  weights <- np:::.np_kernel_weights_direct(
    bws = bws, txdat = training, exdat = evaluation,
    bandwidth.divide = TRUE, kernel.pow = 1.0,
    int.do.tree = np:::DO_TREE_NO
  )
  design <- np:::W.lp(
    xdat = training[, bws$icon, drop = FALSE], degree = degree,
    basis = basis, bernstein.basis = bernstein
  )
  evaluation_design <- np:::W.lp(
    xdat = training[, bws$icon, drop = FALSE],
    exdat = evaluation[, bws$icon, drop = FALSE], degree = degree,
    gradient.vec = derivative, basis = basis,
    bernstein.basis = bernstein
  )

  result <- matrix(NA_real_, nrow(evaluation), nrow(training))
  for (row in seq_len(nrow(evaluation))) {
    weight <- weights[, row]
    coefficient <- solve(
      crossprod(design, design * weight), evaluation_design[row, ]
    )
    result[row, ] <- weight * drop(design %*% coefficient)
  }
  result
}

test_that("beta LP hats agree with independent WLS across orders and bandwidth topologies", {
  set.seed(2026080111L)
  n <- 67L
  training <- data.frame(
    x1 = runif(n, 0.04, 0.96),
    x2 = runif(n, 0.05, 0.95)
  )
  response <- sin(2 * pi * training$x1) + 0.35 * training$x2
  evaluation <- training[c(3L, 19L, 41L, 63L), , drop = FALSE]

  for (order in c(2L, 4L, 6L, 8L)) {
    for (bwtype in c("fixed", "generalized_nn", "adaptive_nn")) {
      bandwidth <- if (identical(bwtype, "fixed")) c(0.32, 0.35) else c(19, 21)
      bernstein <- order %in% c(4L, 8L)
      bws <- npregbw(
        xdat = training, ydat = response, bws = bandwidth,
        bandwidth.compute = FALSE, bwscaling = FALSE,
        regtype = "lp", degree = c(2L, 1L), basis = "glp",
        bernstein.basis = bernstein, bwtype = bwtype,
        ckertype = "beta", ckerorder = order,
        ckerbound = "fixed", ckerlb = c(0, 0), ckerub = c(1, 1)
      )
      actual <- npreghat(
        bws = bws, txdat = training, exdat = evaluation,
        output = "matrix"
      )
      oracle <- beta_lp_hat_oracle(
        bws, training, evaluation, c(2L, 1L),
        bernstein = bernstein
      )
      fit <- npreg(
        bws = bws, txdat = training, tydat = response,
        exdat = evaluation
      )

      expect_equal(
        matrix(as.double(actual), nrow = nrow(actual)), oracle,
        tolerance = 2e-8,
        info = paste("order", order, "bandwidth", bwtype)
      )
      expect_equal(
        drop(actual %*% response), fitted(fit), tolerance = 2e-8,
        info = paste("order", order, "bandwidth", bwtype)
      )
    }
  }
})

test_that("beta LL and raw degree-one LP hats use one canonical route", {
  set.seed(2026080112L)
  training <- data.frame(x1 = runif(73L), x2 = runif(73L))
  response <- cos(training$x1) - training$x2
  common <- list(
    xdat = training, ydat = response, bws = c(0.24, 0.28),
    bandwidth.compute = FALSE, bwscaling = FALSE,
    ckertype = "beta", ckerorder = 6L,
    ckerbound = "fixed", ckerlb = c(0, 0), ckerub = c(1, 1)
  )
  ll <- do.call(npregbw, c(common, list(regtype = "ll")))
  lp <- do.call(npregbw, c(
    common,
    list(regtype = "lp", degree = c(1L, 1L), basis = "glp",
         bernstein.basis = FALSE)
  ))

  expect_identical(
    as.double(npreghat(bws = ll, txdat = training)),
    as.double(npreghat(bws = lp, txdat = training))
  )
})

test_that("mixed beta hats preserve categorical compression exactly", {
  old <- options(np.categorical.compress = FALSE)
  on.exit(options(old), add = TRUE)
  set.seed(2026080113L)
  n <- 72L
  training <- data.frame(
    x = runif(n, 0.03, 0.97),
    u = factor(sample(letters[1:3], n, replace = TRUE)),
    o = ordered(sample(1:4, n, replace = TRUE))
  )
  response <- sin(4 * training$x) + 0.12 * as.integer(training$u) +
    0.07 * as.integer(training$o)
  evaluation <- training[c(2L, 9L, 20L, 43L, 67L), , drop = FALSE]
  bws <- npregbw(
    xdat = training, ydat = response, bws = c(0.18, 0.2, 0.25),
    bandwidth.compute = FALSE, regtype = "lp", degree = 2L,
    basis = "glp", bernstein.basis = FALSE,
    ckertype = "beta", ckerorder = 6L,
    ckerbound = "fixed", ckerlb = 0, ckerub = 1
  )

  options(np.categorical.compress = FALSE)
  dense <- npreghat(bws = bws, txdat = training, exdat = evaluation)
  options(np.categorical.compress = TRUE)
  compressed <- npreghat(bws = bws, txdat = training, exdat = evaluation)
  fit <- npreg(
    bws = bws, txdat = training, tydat = response, exdat = evaluation
  )

  expect_identical(as.double(compressed), as.double(dense))
  expect_equal(drop(compressed %*% response), fitted(fit), tolerance = 2e-9)
})

test_that("beta matrix-free apply shares scalar and general hat arithmetic", {
  old <- options(np.npreghat.apply.memory.threshold.mb = 0)
  on.exit(options(old), add = TRUE)
  set.seed(2026080114L)
  n <- 79L
  training <- data.frame(x1 = runif(n, 0.04, 0.96), x2 = runif(n, 0.04, 0.96))
  response <- sin(training$x1) + training$x2
  rhs <- cbind(response, response^2, cos(response))
  evaluation <- training[c(4L, 18L, 39L, 58L, 74L), , drop = FALSE]

  for (spec in list(
    list(regtype = "lc", degree = c(0L, 0L), bernstein = FALSE),
    list(regtype = "lp", degree = c(2L, 1L), bernstein = TRUE)
  )) {
    bws <- npregbw(
      xdat = training, ydat = response, bws = c(0.31, 0.34),
      bandwidth.compute = FALSE, bwscaling = FALSE,
      regtype = spec$regtype, degree = spec$degree, basis = "glp",
      bernstein.basis = spec$bernstein,
      ckertype = "beta", ckerorder = 8L,
      ckerbound = "fixed", ckerlb = c(0, 0), ckerub = c(1, 1)
    )
    hat <- npreghat(bws = bws, txdat = training, exdat = evaluation)
    applied <- npreghat(
      bws = bws, txdat = training, exdat = evaluation,
      y = rhs, output = "apply"
    )

    expect_equal(
      matrix(as.double(applied), nrow = nrow(applied)),
      matrix(as.double(hat %*% rhs), nrow = nrow(hat)),
      tolerance = 3e-9
    )
  }
})

test_that("beta LP derivative hats and apply share independent WLS arithmetic", {
  old <- options(np.npreghat.apply.memory.threshold.mb = 0)
  on.exit(options(old), add = TRUE)
  set.seed(2026080115L)
  n <- 71L
  training <- data.frame(
    x1 = runif(n, 0.04, 0.96),
    x2 = runif(n, 0.04, 0.96)
  )
  evaluation <- training[c(2L, 13L, 29L, 52L, 68L), , drop = FALSE]
  response <- sin(2 * pi * training$x1) + training$x2^2
  rhs <- cbind(response, response^2, cos(response))

  for (spec in list(
    list(type = "fixed", bandwidth = c(0.31, 0.34), order = 4L,
         degree = c(1L, 2L), bernstein = FALSE),
    list(type = "generalized_nn", bandwidth = c(19, 21), order = 6L,
         degree = c(2L, 1L), bernstein = TRUE),
    list(type = "adaptive_nn", bandwidth = c(19, 21), order = 8L,
         degree = c(2L, 1L), bernstein = FALSE)
  )) {
    bws <- npregbw(
      xdat = training, ydat = response, bws = spec$bandwidth,
      bandwidth.compute = FALSE, bwscaling = FALSE,
      regtype = "lp", degree = spec$degree, basis = "glp",
      bernstein.basis = spec$bernstein, bwtype = spec$type,
      ckertype = "beta", ckerorder = spec$order,
      ckerbound = "fixed", ckerlb = c(0, 0), ckerub = c(1, 1)
    )
    actual <- npreghat(
      bws = bws, txdat = training, exdat = evaluation,
      s = c(1L, 0L), output = "matrix"
    )
    oracle <- beta_lp_hat_oracle(
      bws, training, evaluation, spec$degree,
      bernstein = spec$bernstein, derivative = c(1L, 0L)
    )
    applied <- npreghat(
      bws = bws, txdat = training, exdat = evaluation,
      y = rhs, s = c(1L, 0L), output = "apply"
    )

    expect_equal(
      matrix(as.double(actual), nrow = nrow(actual)), oracle,
      tolerance = 3e-8
    )
    expect_equal(
      matrix(as.double(applied), nrow = nrow(applied)),
      matrix(as.double(actual %*% rhs), nrow = nrow(actual)),
      tolerance = 3e-8
    )
  }
})

test_that("beta matrix-free apply flushes a bounded partial influence block", {
  old <- options(np.npreghat.apply.memory.threshold.mb = 0)
  on.exit(options(old), add = TRUE)
  set.seed(2026080116L)
  n <- 269L
  training <- data.frame(x = runif(n, 0.03, 0.97))
  response <- sin(5 * training$x)
  rhs <- cbind(response, response^2, cos(response))
  bws <- npregbw(
    xdat = training, ydat = response, bws = 0.22,
    bandwidth.compute = FALSE, bwscaling = FALSE,
    regtype = "lp", degree = 2L, basis = "glp",
    bernstein.basis = TRUE, ckertype = "beta", ckerorder = 6L,
    ckerbound = "fixed", ckerlb = 0, ckerub = 1
  )
  hat <- npreghat(bws = bws, txdat = training, output = "matrix")
  applied <- npreghat(
    bws = bws, txdat = training, y = rhs, output = "apply"
  )

  expect_equal(
    matrix(as.double(applied), nrow = nrow(applied)),
    matrix(as.double(hat %*% rhs), nrow = nrow(hat)),
    tolerance = 3e-9
  )
})

test_that("common-scaled beta hats survive complete absolute-weight underflow", {
  old <- options(np.npreghat.apply.memory.threshold.mb = 0)
  on.exit(options(old), add = TRUE)
  set.seed(6219)
  n <- 100L
  p <- 32L
  log_complement <- matrix(runif(n * p, -1.4, -0.8), nrow = n)
  log_complement <- log_complement - rowMeans(log_complement) - 1.1
  training <- as.data.frame(1.0 - exp(log_complement))
  names(training) <- paste0("x", seq_len(p))
  coefficient <- seq_len(p) / p
  response <- 1.0 + drop(as.matrix(training) %*% coefficient)
  evaluation <- as.data.frame(matrix(0.0, nrow = 1L, ncol = p))
  names(evaluation) <- names(training)

  for (regtype in c("lc", "ll")) {
    bws <- npregbw(
      xdat = training, ydat = response, bws = rep.int(0.1, p),
      bandwidth.compute = FALSE, bwscaling = FALSE, regtype = regtype,
      ckertype = "beta", ckerorder = 8L,
      ckerbound = "fixed", ckerlb = rep.int(0.0, p),
      ckerub = rep.int(1.0, p)
    )
    absolute <- npksum(
      bws = bws, txdat = training, exdat = evaluation,
      return.kernel.weights = TRUE
    )
    hat <- npreghat(bws = bws, txdat = training, exdat = evaluation)
    fit <- npreg(
      bws = bws, txdat = training, tydat = response,
      exdat = evaluation
    )
    applied <- npreghat(
      bws = bws, txdat = training, exdat = evaluation,
      y = cbind(response, response^2), output = "apply"
    )

    expect_true(all(absolute$kw == 0.0))
    expect_true(all(is.finite(hat)))
    expect_equal(rowSums(hat), 1.0, tolerance = 2e-8)
    expect_equal(drop(hat %*% response), fitted(fit), tolerance = 2e-8)
    expect_equal(
      matrix(as.double(applied), nrow = nrow(applied)),
      matrix(
        as.double(hat %*% cbind(response, response^2)),
        nrow = nrow(hat)
      ),
      tolerance = 2e-8
    )
  }
})

test_that("beta LP apply retains bounded row storage and no obsolete owner", {
  source <- paste(
    readLines(npRmpi_test_source_path("src", "jksum.c"), warn = FALSE),
    collapse = "\n"
  )
  start <- regexpr("int np_regression_lp_apply_matrix(", source, fixed = TRUE)
  stop <- regexpr(
    "static void np_conditional_yrow_ctx_clear(", source, fixed = TRUE
  )
  expect_gt(start, 0L)
  expect_gt(stop, start)
  apply_source <- substr(source, start, stop - 1L)

  expect_match(apply_source, "8U * 1024U * 1024U", fixed = TRUE)
  expect_match(
    apply_source, "np_reghat_lp_workspace_influence_row", fixed = TRUE
  )
  expect_false(grepl(
    "np_beta_regression_lp_outer_row_canonical", apply_source, fixed = TRUE
  ))
  expect_false(grepl(
    "num_eval*(size_t)num_train*sizeof(double)", apply_source, fixed = TRUE
  ))
})
