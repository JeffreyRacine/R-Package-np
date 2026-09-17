h6_extract_lp_owner <- function(source) {
  start <- gregexpr(
    "static SEXP np_regression_general_lp_fit_execute(void *data)",
    source, fixed = TRUE
  )[[1L]]
  finish <- gregexpr(
    "static int np_regression_general_lp_fit(\n",
    source, fixed = TRUE
  )[[1L]]
  if (length(start) != 1L || start[[1L]] < 1L ||
      length(finish) != 1L || finish[[1L]] < 1L)
    stop("H6 owner requires unique start and finish anchors", call. = FALSE)
  if (finish[[1L]] <= start[[1L]])
    stop("H6 owner anchors are out of order", call. = FALSE)
  substr(source, start[[1L]], finish[[1L]] - 1L)
}

h6_explicit_lp_bw <- function(xdat, ydat, bws, bwtype = "fixed",
                              degree = 2L, basis = "glp",
                              bernstein = FALSE, ...) {
  npregbw(
    xdat = xdat,
    ydat = ydat,
    bws = bws,
    bandwidth.compute = FALSE,
    bwscaling = FALSE,
    bwtype = bwtype,
    regtype = "lp",
    degree = degree,
    degree.select = "manual",
    basis = basis,
    bernstein.basis = bernstein,
    ...
  )
}

h6_hc0_lp_oracle <- function(bws, txdat, tydat, exdat = NULL,
                              gradient.order = NULL) {
  response <- as.double(tydat)
  training.hat <- unclass(suppressWarnings(npreghat(
    bws = bws, txdat = txdat, output = "matrix"
  )))
  evaluation.hat <- if (is.null(exdat)) {
    training.hat
  } else {
    unclass(suppressWarnings(npreghat(
      bws = bws, txdat = txdat, exdat = exdat, output = "matrix"
    )))
  }
  residual <- hc0_normalized_training_residual(training.hat, response)
  if (is.null(gradient.order))
    gradient.order <- rep.int(1L, bws$ncon)

  n.eval <- nrow(evaluation.hat)
  gradient <- matrix(NA_real_, nrow = n.eval, ncol = bws$ncon)
  gradient.stderr <- matrix(NA_real_, nrow = n.eval, ncol = bws$ncon)
  derivative.hat <- vector("list", bws$ncon)
  for (coordinate in seq_len(bws$ncon)) {
    derivative <- integer(bws$ncon)
    derivative[[coordinate]] <- gradient.order[[coordinate]]
    args <- list(
      bws = bws,
      txdat = txdat,
      s = derivative,
      output = "matrix"
    )
    if (!is.null(exdat))
      args$exdat <- exdat
    derivative.hat[[coordinate]] <-
      unclass(suppressWarnings(do.call(npreghat, args)))
    gradient[, coordinate] <-
      drop(derivative.hat[[coordinate]] %*% response)
    gradient.stderr[, coordinate] <- sqrt(drop(
      (derivative.hat[[coordinate]]^2) %*% (residual^2)
    ))
  }

  list(
    mean.stderr = sqrt(drop((evaluation.hat^2) %*% (residual^2))),
    gradient = gradient,
    gradient.stderr = gradient.stderr,
    residual = residual,
    derivative.hat = derivative.hat
  )
}

h6_actual_lp_oracle <- function(bws, txdat, tydat, exdat = NULL,
                                gradient.order = NULL,
                                unit.response = FALSE,
                                constant.reproduction = FALSE) {
  response <- as.double(tydat)
  if (is.null(gradient.order))
    gradient.order <- rep.int(1L, bws$ncon)
  n.train <- nrow(txdat)
  n.eval <- if (is.null(exdat)) n.train else nrow(exdat)
  mean.hat <- matrix(NA_real_, nrow = n.eval, ncol = n.train)
  training.hat <- matrix(NA_real_, nrow = n.train, ncol = n.train)
  derivative.hat <- lapply(
    seq_len(bws$ncon),
    function(...) matrix(NA_real_, nrow = n.eval, ncol = n.train)
  )

  training.fit <- npreg(
    bws = bws, txdat = txdat, tydat = response,
    gradients = FALSE, se = FALSE
  )
  evaluation.args <- list(
    bws = bws,
    txdat = txdat,
    tydat = response,
    gradients = TRUE,
    gradient.order = gradient.order,
    se = FALSE
  )
  if (!is.null(exdat))
    evaluation.args$exdat <- exdat
  evaluation.fit <- suppressWarnings(do.call(npreg, evaluation.args))
  for (donor in seq_len(n.train)) {
    if (isTRUE(unit.response)) {
      # Avoid subtracting nearly equal fits before normalizing a near-identity
      # residual map. This still exercises the public forward solve at the
      # fixed bandwidth/degree; it does not substitute a different hat owner.
      unit <- numeric(n.train)
      unit[[donor]] <- 1
      unit.args <- evaluation.args
      unit.args$tydat <- unit
      unit.fit <- suppressWarnings(do.call(npreg, unit.args))
      mean.hat[, donor] <- unit.fit$mean
      training.hat[, donor] <- if (is.null(exdat)) {
        unit.fit$mean
      } else {
        npreg(bws = bws, txdat = txdat, tydat = unit,
              gradients = FALSE, se = FALSE)$mean
      }
      for (coordinate in seq_len(bws$ncon))
        derivative.hat[[coordinate]][, donor] <-
          unit.fit$grad[, which(bws$icon)[[coordinate]]]
      next
    }
    plus <- response
    minus <- response
    plus[[donor]] <- plus[[donor]] + 1
    minus[[donor]] <- minus[[donor]] - 1
    plus.args <- evaluation.args
    minus.args <- evaluation.args
    plus.args$tydat <- plus
    minus.args$tydat <- minus
    plus.fit <- suppressWarnings(do.call(npreg, plus.args))
    minus.fit <- suppressWarnings(do.call(npreg, minus.args))
    mean.hat[, donor] <- (plus.fit$mean - minus.fit$mean) / 2
    if (is.null(exdat)) {
      training.hat[, donor] <- mean.hat[, donor]
    } else {
      # Retain the actual public ridge map, not a substituted hat owner.
      plus.training <- npreg(
        bws = bws, txdat = txdat, tydat = plus,
        gradients = FALSE, se = FALSE
      )
      minus.training <- npreg(
        bws = bws, txdat = txdat, tydat = minus,
        gradients = FALSE, se = FALSE
      )
      training.hat[, donor] <- (plus.training$mean - minus.training$mean) / 2
    }
    for (coordinate in seq_len(bws$ncon))
      derivative.hat[[coordinate]][, donor] <-
        (plus.fit$grad[, which(bws$icon)[[coordinate]]] -
           minus.fit$grad[, which(bws$icon)[[coordinate]]]) / 2
  }

  residual <- hc0_normalized_training_residual(
    training.hat, response, training.mean = training.fit$mean,
    constant.reproduction = constant.reproduction
  )
  gradient <- evaluation.fit$grad[, which(bws$icon), drop = FALSE]
  gradient.stderr <- vapply(
    derivative.hat,
    function(hat) sqrt(drop((hat^2) %*% (residual^2))),
    numeric(n.eval)
  )
  if (bws$ncon == 1L)
    gradient.stderr <- matrix(gradient.stderr, ncol = 1L)

  list(
    mean.stderr = sqrt(drop((mean.hat^2) %*% (residual^2))),
    gradient = gradient,
    gradient.stderr = gradient.stderr,
    residual = residual,
    derivative.hat = derivative.hat
  )
}

h6_expect_lp_derivative_contract <- function(
    bws, txdat, tydat, exdat = NULL, gradient.order = NULL,
    tolerance = 3e-9, oracle.method = c("hat", "actual")) {
  oracle.method <- match.arg(oracle.method)
  args <- list(
    bws = bws,
    txdat = txdat,
    tydat = tydat,
    gradients = TRUE
  )
  if (!is.null(exdat))
    args$exdat <- exdat
  if (!is.null(gradient.order))
    args$gradient.order <- gradient.order

  without.se <- do.call(npreg, c(args, list(se = FALSE)))
  with.se <- do.call(npreg, c(args, list(se = TRUE)))
  oracle <- if (identical(oracle.method, "hat")) {
    h6_hc0_lp_oracle(
      bws, txdat, tydat, exdat, gradient.order = gradient.order
    )
  } else {
    h6_actual_lp_oracle(
      bws, txdat, tydat, exdat, gradient.order = gradient.order,
      # Beta LP uses the accepted constant-completed residual map. The
      # separate literal actual-ridge comparator remains in place elsewhere.
      unit.response = identical(bws$ckertype, "beta"),
      constant.reproduction = identical(bws$ckertype, "beta")
    )
  }

  expect_identical(with.se$mean, without.se$mean)
  expect_identical(with.se$grad, without.se$grad)
  expect_identical(with.se$xtra, without.se$xtra)
  expect_equal(with.se$merr, oracle$mean.stderr, tolerance = tolerance)
  expect_equal(
    with.se$grad[, bws$icon, drop = FALSE],
    oracle$gradient,
    tolerance = tolerance
  )
  expect_equal(
    with.se$gerr[, bws$icon, drop = FALSE],
    oracle$gradient.stderr,
    tolerance = tolerance
  )
  expect_true(all(is.finite(with.se$gerr[, bws$icon, drop = FALSE])))
  expect_true(all(with.se$gerr[, bws$icon, drop = FALSE] >= 0))

  categorical <- setdiff(seq_len(ncol(with.se$gerr)), which(bws$icon))
  if (length(categorical)) {
    expect_true(all(is.finite(with.se$gerr[, categorical, drop = FALSE])))
    expect_true(all(with.se$gerr[, categorical, drop = FALSE] >= 0))
  }
  invisible(list(fit = with.se, oracle = oracle))
}

test_that("general-LP continuous derivative HC0 matches actual hats across bandwidth modes", {
  old <- options(np.messages = FALSE, np.tree = TRUE)
  on.exit(options(old), add = TRUE)

  n <- 31L
  txdat <- data.frame(
    x1 = seq(-1.1, 1.15, length.out = n),
    x2 = 0.55 * cos(seq(-0.6, 2.4, length.out = n)),
    u = factor(rep(c("a", "b", "c"), length.out = n))
  )
  tydat <- sin(1.3 * txdat$x1) + 0.45 * txdat$x2 +
    c(a = -0.2, b = 0.15, c = 0.4)[txdat$u] + seq_len(n) / 130
  exdat <- txdat[c(2L, 8L, 15L, 23L, 30L), , drop = FALSE]

  for (bwtype in c("fixed", "generalized_nn", "adaptive_nn")) {
    bandwidth <- if (identical(bwtype, "fixed")) {
      c(0.46, 0.39, 0.2)
    } else {
      c(9, 9, 0.2)
    }
    bw <- h6_explicit_lp_bw(
      txdat, tydat, bandwidth, bwtype = bwtype,
      degree = c(2L, 1L), basis = "glp", bernstein = FALSE
    )
    h6_expect_lp_derivative_contract(bw, txdat, tydat)
    h6_expect_lp_derivative_contract(bw, txdat, tydat, exdat)
  }
})

test_that("all LP bases and requested derivative orders share the HC0 adjoint block", {
  old <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old), add = TRUE)

  n <- 29L
  txdat <- data.frame(
    x1 = seq(0.03, 0.97, length.out = n),
    x2 = (((seq_len(n) * 7L) %% n) + 0.5) / n
  )
  tydat <- sin(2 * pi * txdat$x1) + 0.3 * txdat$x2^2 +
    seq_len(n) / 150
  exdat <- txdat[c(3L, 10L, 18L, 27L), , drop = FALSE]
  specifications <- list(
    list(basis = "glp", bernstein = FALSE),
    list(basis = "glp", bernstein = TRUE),
    list(basis = "tensor", bernstein = TRUE)
  )

  for (specification in specifications) {
    bw <- h6_explicit_lp_bw(
      txdat, tydat, c(0.3, 0.33), degree = c(2L, 2L),
      basis = specification$basis,
      bernstein = specification$bernstein,
      ckerbound = "fixed", ckerlb = c(0, 0), ckerub = c(1, 1)
    )
    h6_expect_lp_derivative_contract(
      bw, txdat, tydat, exdat, gradient.order = c(2L, 1L),
      tolerance = 2e-8, oracle.method = "actual"
    )
  }
})

test_that("beta general-LP derivative HC0 covers every order and bandwidth mode", {
  old <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old), add = TRUE)

  n <- 15L
  txdat <- data.frame(
    x1 = seq(0.04, 0.96, length.out = n),
    x2 = 0.5 + 0.43 * sin(seq(0.2, 2.8, length.out = n))
  )
  tydat <- sin(4 * txdat$x1) + 0.35 * txdat$x2 + seq_len(n) / 140
  exdat <- txdat[c(2L, 5L, 8L, 12L, 15L), , drop = FALSE]

  for (order in c(2L, 4L, 6L, 8L)) {
    for (bwtype in c("fixed", "generalized_nn", "adaptive_nn")) {
      cell.x <- txdat
      cell.y <- tydat
      cell.ex <- exdat
      if (order == 2L && identical(bwtype, "adaptive_nn")) {
        # The original curved design is nearly interpolating at its last
        # row (q about 2.7e-13). Public unit responses change beta's retained
        # moment scale/forward rounding before q normalization, so that map
        # is not a portable strict numerical oracle for the original stress
        # fixture below. Retain the same kernel/order/bandwidth/basis contract
        # on a non-nearly-interpolating design for this reference comparison.
        cell.x$x2 <- (((seq_len(n) * 7L) %% n) + 0.5) / n
        cell.y <- sin(4 * cell.x$x1) + 0.35 * cell.x$x2 + seq_len(n) / 140
        cell.ex <- cell.x[c(2L, 5L, 8L, 12L, 15L), , drop = FALSE]
      }
      bw <- h6_explicit_lp_bw(
        cell.x, cell.y,
        if (identical(bwtype, "fixed")) c(0.27, 0.3) else c(7, 7),
        bwtype = bwtype, degree = c(2L, 1L),
        basis = "glp", bernstein = order %in% c(4L, 8L),
        ckertype = "beta", ckerorder = order,
        ckerbound = "fixed", ckerlb = c(0, 0), ckerub = c(1, 1)
      )
      h6_expect_lp_derivative_contract(
        bw, cell.x, cell.y, cell.ex, tolerance = 1e-5,
        oracle.method = "actual"
      )
    }
  }
})

test_that("near-interpolating beta ANN retains finite requested inference", {
  old <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old), add = TRUE)
  n <- 15L
  x <- data.frame(x1 = seq(0.04, 0.96, length.out = n),
    x2 = 0.5 + 0.43 * sin(seq(0.2, 2.8, length.out = n)))
  y <- sin(4 * x$x1) + 0.35 * x$x2 + seq_len(n) / 140
  ex <- x[c(2L, 5L, 8L, 12L, 15L), , drop = FALSE]
  bw <- h6_explicit_lp_bw(x, y, c(7, 7), bwtype = "adaptive_nn",
    degree = c(2L, 1L), basis = "glp", bernstein = FALSE,
    ckertype = "beta", ckerorder = 2L, ckerbound = "fixed",
    ckerlb = c(0, 0), ckerub = c(1, 1))
  point <- npreg(bws = bw, txdat = x, tydat = y, exdat = ex,
                 gradients = TRUE, se = FALSE)
  inference <- npreg(bws = bw, txdat = x, tydat = y, exdat = ex,
                     gradients = TRUE, se = TRUE)
  expect_identical(inference$mean, point$mean)
  expect_identical(inference$grad, point$grad)
  expect_identical(inference$xtra, point$xtra)
  expect_true(all(is.finite(inference$merr)))
  expect_true(all(is.finite(inference$gerr)))
  expect_true(all(inference$merr >= 0))
  expect_true(all(inference$gerr >= 0))
  # Exact retained adjoint/basis/kernel inputs for this original fixture are
  # separately checked by the campaign's 90-digit donor-covariance oracle.
  # Do not replace that proof with a hard-coded platform-dependent SE or a
  # relaxed tolerance against a response-perturbed public forward map.
})

test_that("general-LP derivative HC0 follows accepted ridge and all-large maps", {
  old <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old), add = TRUE)

  ridge.x <- data.frame(
    x = c(-1, -1, -0.35, -0.35, 0.2, 0.2, 0.75, 0.75, 1.1)
  )
  ridge.y <- c(0.2, 0.35, -0.1, 0.15, 0.8, 0.65, 1.4, 1.1, 1.7)
  ridge.eval <- data.frame(x = c(-0.9, -0.25, 0.3, 0.9))
  ridge.bw <- suppressWarnings(h6_explicit_lp_bw(
    ridge.x, ridge.y, 0.19, degree = 3L,
    ckertype = "epanechnikov"
  ))
  h6_expect_lp_derivative_contract(
    ridge.bw, ridge.x, ridge.y, ridge.eval, tolerance = 3e-8,
    oracle.method = "actual"
  )

  large.x <- data.frame(x = seq(-1, 1, length.out = 21L))
  large.y <- sin(large.x$x) + seq_len(nrow(large.x)) / 90
  large.eval <- data.frame(x = c(-0.8, -0.15, 0.45, 0.9))
  large.bw <- h6_explicit_lp_bw(
    large.x, large.y, 1e16, degree = 2L
  )
  h6_expect_lp_derivative_contract(
    large.bw, large.x, large.y, large.eval, tolerance = 2e-9,
    oracle.method = "actual"
  )
})

test_that("balanced q-adjusted multiplier covariance equals general-LP HC0", {
  old <- options(np.messages = FALSE, np.tree = FALSE)
  on.exit(options(old), add = TRUE)

  n <- 16L
  txdat <- data.frame(x = seq(-1, 1, length.out = n))
  tydat <- sin(1.7 * txdat$x) + seq_len(n) / 60
  exdat <- data.frame(x = c(-0.75, -0.2, 0.35, 0.8))
  bw <- h6_explicit_lp_bw(txdat, tydat, 0.48, degree = 2L)
  result <- h6_expect_lp_derivative_contract(bw, txdat, tydat, exdat)

  hadamard <- matrix(1, 1L, 1L)
  while (nrow(hadamard) < n)
    hadamard <- rbind(
      cbind(hadamard, hadamard),
      cbind(hadamard, -hadamard)
    )
  multipliers <- rbind(hadamard, -hadamard)
  # This is an algebraic covariance identity using q-adjusted donors, not a
  # claim that the separate public wild-bootstrap residual law was changed.
  derivative.donor <- sweep(
    result$oracle$derivative.hat[[1L]],
    2L,
    result$oracle$residual,
    `*`
  )
  bootstrap.deviation <- derivative.donor %*% t(multipliers)
  bootstrap.stderr <- sqrt(rowMeans(bootstrap.deviation^2))

  expect_equal(
    result$fit$gerr[, bw$icon], bootstrap.stderr,
    tolerance = 2e-12
  )
})

test_that("H6 batches adjoint directions without another covariance owner", {
  source <- paste(
    readLines(test_path("..", "..", "src", "jksum.c"), warn = FALSE),
    collapse = "\n"
  )
  owner <- h6_extract_lp_owner(source)

  expect_match(
    owner,
    "if(call->do_grad && call->do_gerr)",
    fixed = TRUE
  )
  expect_match(
    owner,
    "np_lp_solve_workspace_solve_adjoint_factored(",
    fixed = TRUE
  )
  expect_match(
    owner,
    "np_regression_hc0_lp_standard_error_reuse(",
    fixed = TRUE
  )
  # Exact-input reuse is only a computational wrapper. No-cache and miss
  # paths must still use the same donor-square calculation and propagate its
  # failures before publishing a reusable result.
  source.lines <- strsplit(source, "\n", fixed = TRUE)[[1L]]
  reuse <- npRmpi_test_extract_c_function(
    source.lines, "np_regression_hc0_lp_standard_error_reuse")
  donor <- gregexpr("np_regression_hc0_lp_standard_error_with_information(",
                   reuse, fixed = TRUE)[[1L]]
  expect_length(donor[donor > 0L], 2L)
  compact.reuse <- gsub("[[:space:]]+", " ", reuse)
  expect_match(compact.reuse,
    paste0("if(reuse->entry == NULL) return ",
           "np_regression_hc0_lp_standard_error_with_information("),
    fixed = TRUE)
  expect_match(compact.reuse,
    "standard_error, &missing)) return 0; np_inference_reuse_put(",
    fixed = TRUE)
  expect_match(reuse, "if(!np_inference_reuse_get(", fixed = TRUE)
  expect_match(reuse, "if(unavailable != NULL && missing) *unavailable = 1;",
               fixed = TRUE)
  expect_lt(regexpr("np_inference_reuse_get(", reuse, fixed = TRUE)[[1L]],
            donor[[2L]])
  expect_lt(donor[[2L]],
            regexpr("np_inference_reuse_put(", reuse, fixed = TRUE)[[1L]])
  expect_match(gsub("[[:space:]]+", " ", owner),
    paste0("if(ordinary_hc0 && call->do_merr && call->all_large_kernel_rows) ",
           "(void)np_inference_reuse_reserve("), fixed = TRUE)
  cleanup <- npRmpi_test_extract_c_function(
    source.lines, "np_regression_general_lp_fit_owner_cleanup")
  expect_match(cleanup, "np_inference_reuse_clear(&owner->inference_reuse);",
               fixed = TRUE)
  expect_false(grepl(
    "!ordinary_hc0 && call->do_grad && call->do_gerr",
    owner,
    fixed = TRUE
  ))
  expect_false(grepl("npreghat", owner, fixed = TRUE))
  expect_false(grepl("inverse", owner, fixed = TRUE))
})
