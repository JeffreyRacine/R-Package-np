legacy_generalized_lc_derivative_hat <- function(bws, txdat, exdat = NULL, s) {
  no.ex <- is.null(exdat)
  eval.data <- if (no.ex) txdat else exdat
  ntrain <- nrow(txdat)
  neval <- nrow(eval.data)
  target.cont <- which(s == 1L)
  beta.kernel <- identical(bws[["ckertype", exact = TRUE]], "beta")
  ones <- rep.int(1.0, ntrain)
  W <- cbind(diag(ntrain), ones)

  call <- quote(np:::npksum.default(
    bws = bws,
    txdat = txdat,
    exdat = if (no.ex) txdat else eval.data,
    tydat = ones,
    weights = W,
    bandwidth.divide = !beta.kernel,
    permutation.operator = "derivative"
  ))
  out <- if (beta.kernel) suppressWarnings(eval(call)) else eval(call)
  ks <- out$ksum
  ps <- out$p.ksum
  if (!is.matrix(ks))
    ks <- matrix(ks, nrow = ntrain + 1L, ncol = neval)
  if (length(dim(ps)) == 3L)
    ps <- ps[, , target.cont, drop = TRUE]
  if (!is.matrix(ps))
    ps <- matrix(ps, nrow = ntrain + 1L, ncol = neval)

  sk <- ks[ntrain + 1L, ]
  dsk <- ps[ntrain + 1L, ]
  H <- t(
    ps[seq_len(ntrain), , drop = FALSE] / rep(sk, each = ntrain) -
      ks[seq_len(ntrain), , drop = FALSE] *
        rep(dsk / (sk^2), each = ntrain)
  )
  if (beta.kernel) {
    undefined <- !is.finite(sk) | !is.finite(dsk) |
      apply(!is.finite(ps), 2L, any)
    if (any(undefined)) {
      H[undefined, ] <- NA_real_
      warning(sprintf(
        "beta derivative hat matrix contains %d undefined endpoint row(s)",
        sum(undefined)
      ), call. = FALSE)
    }
  }
  H
}

capture_generalized_hat_contract <- function(expr) {
  warnings <- character()
  value <- withCallingHandlers(
    expr,
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(value = value, warnings = warnings)
}

test_that("generalized-NN lc derivative native rows preserve the fitted statistic", {
  set.seed(2026072261L)
  n <- 65L
  tx <- data.frame(
    x1 = runif(n, -0.7, 0.9),
    x2 = runif(n, -0.4, 1.1),
    f = factor(sample(letters[1:3], n, replace = TRUE))
  )
  ex <- tx[c(seq_len(9L), n), , drop = FALSE]
  y <- sin(tx$x1) - 0.4 * tx$x2 + as.integer(tx$f) / 5
  cases <- list(
    gaussian = list(kernel = "gaussian", tree = FALSE),
    epan_tree = list(kernel = "epanechnikov", tree = TRUE),
    uniform_control = list(kernel = "uniform", tree = TRUE)
  )

  for (case in cases) {
    withr::local_options(list(np.tree = case$tree))
    bw <- suppressWarnings(npregbw(
      xdat = tx, ydat = y, regtype = "lc", bwtype = "generalized_nn",
      ckertype = case$kernel, ckerorder = 6L,
      bws = c(15L, 19L, 0.35), bandwidth.compute = FALSE
    ))
    for (evaluation in list(NULL, ex)) {
      fast.args <- list(
        bws = bw, txdat = tx, output = "matrix", s = c(0L, 1L)
      )
      if (!is.null(evaluation))
        fast.args$exdat <- evaluation
      fast <- capture_generalized_hat_contract(do.call(npreghat, fast.args))
      legacy <- capture_generalized_hat_contract(
        legacy_generalized_lc_derivative_hat(
          bws = bw, txdat = tx, exdat = evaluation, s = c(0L, 1L)
        )
      )
      fit <- if (is.null(evaluation)) {
        npreg(bws = bw, gradients = TRUE)
      } else {
        npreg(bws = bw, exdat = evaluation, gradients = TRUE)
      }
      expect_equal(drop(fast$value %*% y), fit$grad[, 2L],
                   tolerance = 1e-12)
      if (!is.null(evaluation)) {
        expect_equal(as.double(fast$value), as.double(legacy$value),
                     tolerance = 1e-13)
      }
      expect_identical(fast$warnings, legacy$warnings)
    }
  }
})

test_that("generalized-NN beta derivative route preserves endpoint behavior", {
  tx <- data.frame(x = c(0, 0.01, 0.04, seq(0.08, 0.92, length.out = 35L),
                           0.96, 0.99, 1))
  y <- sin(2 * pi * tx$x)
  bw <- npregbw(
    xdat = tx, ydat = y, bws = 11L,
    bandwidth.compute = FALSE, regtype = "lc", bwtype = "generalized_nn",
    ckertype = "beta", ckerorder = 6L,
    ckerbound = "fixed", ckerlb = 0, ckerub = 1
  )
  fast <- capture_generalized_hat_contract(npreghat(
    bws = bw, txdat = tx, output = "matrix", s = 1L
  ))
  legacy <- capture_generalized_hat_contract(
    legacy_generalized_lc_derivative_hat(bws = bw, txdat = tx, s = 1L)
  )
  # Training rows use self-excluded NN radii; the external-query legacy hat
  # above intentionally has a different radius. Retain its endpoint warnings.
  fit <- suppressWarnings(npreg(bws=bw, txdat=tx, tydat=y, gradients=TRUE))
  interior <- tx$x>0 & tx$x<1
  expect_equal(drop(fast$value %*% y)[interior],
               as.double(fit$grad)[interior], tolerance=1e-12)
  expect_identical(is.na(fast$value), is.na(legacy$value))
  # Independent derivative of normalized beta weights with the NN radius fixed.
  radius <- vapply(seq_len(nrow(tx)), function(i)
    sort(abs(tx$x[-i]-tx$x[i]))[11L], numeric(1L))
  oracle <- vapply(seq_len(nrow(tx)), function(i) {
    target <- tx$x[i]; donor <- tx$x
    if (target==0 || target==1) return(NA_real_)
    weights <- derivatives <- numeric(length(donor))
    interior <- donor>0 & donor<1
    for (j in 1:3) {
      concentration <- 1/(j*radius[i]^2)
      alpha <- 1+target*concentration
      beta <- 1+(1-target)*concentration
      part <- (-1)^(j+1L)*choose(3L,j)*dbeta(donor,alpha,beta)
      weights <- weights+part
      derivatives[interior] <- derivatives[interior]+part[interior]*concentration*
        (log(donor[interior])-log1p(-donor[interior])-digamma(alpha)+digamma(beta))
    }
    (sum(derivatives*y)-sum(weights*y)*sum(derivatives)/sum(weights))/sum(weights)
  }, numeric(1L))
  interior <- is.finite(oracle)
  expect_equal(drop(fast$value %*% y)[interior], oracle[interior], tolerance=1e-10)
  expect_identical(fast$warnings, legacy$warnings)
})

test_that("private divided-row seam is opt-in and leaves public rows raw", {
  set.seed(2026072262L)
  tx <- data.frame(x1 = runif(31L), x2 = runif(31L))
  y <- sin(tx$x1) + tx$x2
  bw <- npregbw(
    xdat = tx, ydat = y, regtype = "lc", bwtype = "fixed",
    bws = c(0.2, 0.3), bandwidth.compute = FALSE
  )
  public <- np:::npksum.default(
    bws = bw, txdat = tx, bandwidth.divide = TRUE,
    return.kernel.weights = TRUE
  )
  divided <- np:::npksum.default(
    bws = bw, txdat = tx, bandwidth.divide = TRUE,
    return.kernel.weights = TRUE,
    .np.internal.bandwidth.divide.weights = TRUE
  )
  expect_identical(public$ksum, divided$ksum)
  expect_identical(public$kw / prod(bw$bw[bw$icon]), divided$kw)
  expect_error(
    np:::npksum.default(
      bws = bw, txdat = tx, bandwidth.divide = TRUE,
      .np.internal.bandwidth.divide.weights = TRUE
    ),
    "invalid use of the internal divided kernel-weight route"
  )
  expect_error(
    np:::npksum.default(
      bws = bw, txdat = tx, return.kernel.weights = TRUE,
      .np.internal.bandwidth.divide.weights = TRUE
    ),
    "invalid use of the internal divided kernel-weight route"
  )
})
