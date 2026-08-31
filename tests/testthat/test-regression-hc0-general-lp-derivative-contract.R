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
  residual <- response - drop(training.hat %*% response)
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
                                gradient.order = NULL) {
  response <- as.double(tydat)
  if (is.null(gradient.order))
    gradient.order <- rep.int(1L, bws$ncon)
  n.train <- nrow(txdat)
  n.eval <- if (is.null(exdat)) n.train else nrow(exdat)
  mean.hat <- matrix(NA_real_, nrow = n.eval, ncol = n.train)
  derivative.hat <- lapply(
    seq_len(bws$ncon),
    function(...) matrix(NA_real_, nrow = n.eval, ncol = n.train)
  )

  training.fit <- npreg(
    bws = bws, txdat = txdat, tydat = response,
    gradients = FALSE, se = FALSE
  )
  residual <- response - training.fit$mean
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
    for (coordinate in seq_len(bws$ncon))
      derivative.hat[[coordinate]][, donor] <-
        (plus.fit$grad[, which(bws$icon)[[coordinate]]] -
           minus.fit$grad[, which(bws$icon)[[coordinate]]]) / 2
  }

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
      bws, txdat, tydat, exdat, gradient.order = gradient.order
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
      bw <- h6_explicit_lp_bw(
        txdat, tydat,
        if (identical(bwtype, "fixed")) c(0.27, 0.3) else c(7, 7),
        bwtype = bwtype, degree = c(2L, 1L),
        basis = "glp", bernstein = order %in% c(4L, 8L),
        ckertype = "beta", ckerorder = order,
        ckerbound = "fixed", ckerlb = c(0, 0), ckerub = c(1, 1)
      )
      h6_expect_lp_derivative_contract(
        bw, txdat, tydat, exdat, tolerance = 1e-5,
        oracle.method = "actual"
      )
    }
  }
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

test_that("balanced fixed-design wild covariance equals general-LP HC0", {
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
  start <- regexpr(
    "static SEXP np_regression_general_lp_fit_execute(void *data)",
    source,
    fixed = TRUE
  )[[1L]]
  finish <- regexpr(
    "static int np_regression_general_lp_fit(\n",
    substring(source, start),
    fixed = TRUE
  )[[1L]]
  expect_gt(start, 0L)
  expect_gt(finish, 0L)
  owner <- substr(source, start, start + finish - 2L)

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
    "np_regression_hc0_lp_standard_error(",
    fixed = TRUE
  )
  expect_false(grepl(
    "!ordinary_hc0 && call->do_grad && call->do_gerr",
    owner,
    fixed = TRUE
  ))
  expect_false(grepl("npreghat", owner, fixed = TRUE))
  expect_false(grepl("inverse", owner, fixed = TRUE))
})
