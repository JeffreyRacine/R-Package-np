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

