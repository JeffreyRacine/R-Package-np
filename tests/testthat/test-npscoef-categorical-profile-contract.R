npscoef_cat_profile_case <- function(kind = c("uno", "ord", "mixed", "rly"),
                                     eval = FALSE,
                                     n = 300L,
                                     seed = 1L) {
  kind <- match.arg(kind)
  set.seed(seed)
  xdat <- data.frame(x = rnorm(n))
  if (identical(kind, "uno")) {
    zdat <- data.frame(z = factor(sample(letters[1:4], n, TRUE)))
  } else if (identical(kind, "ord")) {
    zdat <- data.frame(o = ordered(sample(1:5, n, TRUE)))
  } else if (identical(kind, "mixed")) {
    zdat <- data.frame(
      z = factor(sample(letters[1:3], n, TRUE)),
      o = ordered(sample(1:4, n, TRUE))
    )
  } else {
    zdat <- data.frame(o = ordered(sample(1:5, n, TRUE)))
  }

  ydat <- 1 + 0.7 * xdat$x +
    rowSums(as.data.frame(lapply(zdat, as.numeric))) +
    rnorm(n, sd = 0.15)
  bw <- npscoefbw(
    xdat = xdat,
    ydat = ydat,
    zdat = zdat,
    bws = rep(0.25, ncol(zdat)),
    bandwidth.compute = FALSE,
    regtype = "lc",
    okertype = if (identical(kind, "rly")) "racineliyan" else "liracine"
  )
  args <- list(
    bws = bw,
    txdat = xdat,
    tydat = ydat,
    tzdat = zdat,
    se = TRUE,
    iterate = FALSE,
    betas = TRUE
  )
  if (eval) {
    m <- 41L
    args$exdat <- data.frame(x = seq(min(xdat$x), max(xdat$x), length.out = m))
    args$ezdat <- zdat[seq_len(m), , drop = FALSE]
  }

  options(np.tree = FALSE, np.categorical.compress = FALSE)
  dense <- do.call(npscoef, args)
  options(np.tree = FALSE, np.categorical.compress = TRUE)
  profile <- do.call(npscoef, args)
  list(dense = dense, profile = profile)
}

test_that("npscoef all-categorical profile route preserves fitted and errors", {
  old <- options(np.messages = FALSE, np.tree = FALSE,
                 np.categorical.compress = FALSE)
  on.exit(options(old), add = TRUE)

  cases <- c("uno", "ord", "mixed", "rly")
  for (i in seq_along(cases)) {
    out <- npscoef_cat_profile_case(cases[[i]], eval = FALSE,
                                    seed = 20260516L + i)
    expect_equal(out$profile$mean, out$dense$mean, tolerance = 1e-8,
                 info = cases[[i]])
    expect_equal(out$profile$beta, out$dense$beta, tolerance = 1e-8,
                 info = cases[[i]])
    expect_equal(out$profile$merr, out$dense$merr, tolerance = 1e-8,
                 info = cases[[i]])
    expect_equal(out$profile$gerr, out$dense$gerr, tolerance = 1e-8,
                 info = cases[[i]])
  }
})

test_that("npscoef all-categorical profile route preserves bandwidth CV", {
  old <- options(np.messages = FALSE, np.tree = FALSE,
                 np.categorical.compress = FALSE)
  on.exit(options(old), add = TRUE)

  set.seed(20261017L)
  n <- 100L
  dat <- data.frame(
    y = rnorm(n),
    x1 = rnorm(n),
    x2 = rnorm(n),
    z1 = factor(sample(letters[1:3], n, TRUE)),
    z2 = ordered(sample(1:4, n, TRUE)),
    z3 = factor(sample(c("a", "b"), n, TRUE))
  )
  dat$y <- 1 + 0.6 * dat$x1 - 0.3 * dat$x2 +
    as.numeric(dat$z1) - 0.5 * as.numeric(dat$z2) +
    rnorm(n, sd = 0.2)

  options(np.tree = FALSE, np.categorical.compress = FALSE)
  dense <- npscoefbw(y ~ x1 + x2 | z1 + z2 + z3,
                     data = dat, regtype = "lc", nmulti = 1)
  options(np.tree = FALSE, np.categorical.compress = TRUE)
  profile <- npscoefbw(y ~ x1 + x2 | z1 + z2 + z3,
                       data = dat, regtype = "lc", nmulti = 1)

  expect_equal(profile$fval, dense$fval, tolerance = 1e-8)
  expect_gt(profile$num.feval.fast, 0L)
})

test_that("npscoef all-categorical profile route preserves evaluation output", {
  old <- options(np.messages = FALSE, np.tree = FALSE,
                 np.categorical.compress = FALSE)
  on.exit(options(old), add = TRUE)

  cases <- c("uno", "ord", "mixed", "rly")
  for (i in seq_along(cases)) {
    out <- npscoef_cat_profile_case(cases[[i]], eval = TRUE,
                                    seed = 20260616L + i)
    expect_equal(out$profile$mean, out$dense$mean, tolerance = 1e-8,
                 info = cases[[i]])
    expect_equal(out$profile$beta, out$dense$beta, tolerance = 1e-8,
                 info = cases[[i]])
    expect_equal(out$profile$merr, out$dense$merr, tolerance = 1e-8,
                 info = cases[[i]])
    expect_equal(out$profile$gerr, out$dense$gerr, tolerance = 1e-8,
                 info = cases[[i]])
  }
})

test_that("npscoef categorical compression option leaves mixed z dense", {
  old <- options(np.messages = FALSE, np.tree = FALSE,
                 np.categorical.compress = FALSE)
  on.exit(options(old), add = TRUE)

  set.seed(20260716L)
  n <- 120L
  xdat <- data.frame(x = rnorm(n))
  zdat <- data.frame(z = factor(sample(letters[1:3], n, TRUE)),
                     c = runif(n))
  ydat <- 1 + xdat$x + as.numeric(zdat$z) + sin(zdat$c) + rnorm(n, sd = 0.2)
  bw <- npscoefbw(
    xdat = xdat,
    ydat = ydat,
    zdat = zdat,
    bws = c(0.25, 0.3),
    bandwidth.compute = FALSE,
    regtype = "lc"
  )

  options(np.tree = FALSE, np.categorical.compress = FALSE)
  dense <- npscoef(bws = bw, txdat = xdat, tydat = ydat, tzdat = zdat,
                  se = FALSE, iterate = FALSE)
  options(np.tree = FALSE, np.categorical.compress = TRUE)
  profile <- npscoef(bws = bw, txdat = xdat, tydat = ydat, tzdat = zdat,
                    se = FALSE, iterate = FALSE)
  expect_equal(profile$mean, dense$mean, tolerance = 1e-10)
})

test_that("npscoefhat apply uses categorical profile route without changing output", {
  old <- options(np.messages = FALSE, np.tree = FALSE,
                 np.categorical.compress = FALSE)
  on.exit(options(old), add = TRUE)

  set.seed(20260816L)
  n <- 220L
  xdat <- data.frame(x = rnorm(n))
  zdat <- data.frame(
    z = factor(sample(letters[1:4], n, TRUE)),
    o = ordered(sample(1:5, n, TRUE))
  )
  ydat <- 1 + xdat$x + as.numeric(zdat$z) + as.numeric(zdat$o) +
    rnorm(n, sd = 0.2)
  bw <- npscoefbw(
    xdat = xdat,
    ydat = ydat,
    zdat = zdat,
    bws = c(0.25, 0.3),
    bandwidth.compute = FALSE,
    regtype = "lc"
  )
  rhs <- cbind(ydat, ydat^2)

  options(np.tree = FALSE, np.categorical.compress = FALSE)
  dense <- npscoefhat(bws = bw, txdat = xdat, tzdat = zdat,
                     y = rhs, output = "apply")
  options(np.tree = FALSE, np.categorical.compress = TRUE)
  profile <- npscoefhat(bws = bw, txdat = xdat, tzdat = zdat,
                       y = rhs, output = "apply")
  expect_equal(profile, dense, tolerance = 1e-8)
})

test_that("npscoefhat categorical profile ridge uses the canonical weight scale", {
  old <- options(np.messages = FALSE, np.tree = FALSE,
                 np.categorical.compress = TRUE)
  on.exit(options(old), add = TRUE)

  set.seed(20260917L)
  n <- 48L
  xdat <- data.frame(a = runif(n, -1, 1), b = runif(n, -1, 1))
  zdat <- data.frame(
    u = factor(rep(letters[1:3], length.out = n)),
    v = factor(rep(LETTERS[1:2], length.out = n))
  )
  y <- xdat$a - 0.35 * xdat$b + rnorm(n, sd = 0.1)
  exdat <- xdat[1:7, , drop = FALSE]
  ezdat <- zdat[1:7, , drop = FALSE]
  for (ukertype in c("liracine", "aitchisonaitken")) {
    bw <- npscoefbw(
      xdat = xdat, ydat = y, zdat = zdat, bws = c(0.4, 0.3),
      bandwidth.compute = FALSE, ukertype = ukertype
    )
    args <- list(bws = bw, txdat = xdat, tzdat = zdat,
                 exdat = exdat, ezdat = ezdat, ridge = 0.25)

    H <- do.call(npscoefhat, c(args, list(output = "matrix")))
    applied.profile <- do.call(npscoefhat, c(args, list(y = y, output = "apply")))
    options(np.categorical.compress = FALSE)
    applied.direct <- do.call(npscoefhat, c(args, list(y = y, output = "apply")))
    options(np.categorical.compress = TRUE)

    W <- cbind(1, as.matrix(xdat))
    E <- cbind(1, as.matrix(exdat))
    oracle <- vapply(seq_len(nrow(exdat)), function(i) {
      weights <- rep(1, n)
      for (j in seq_along(zdat)) {
        same <- zdat[[j]] == ezdat[[j]][i]
        lambda <- c(0.4, 0.3)[j]
        nc <- nlevels(zdat[[j]])
        weights <- weights * if (ukertype == "liracine")
          ifelse(same, 1, lambda) / (1 + (nc - 1) * lambda) else
          ifelse(same, 1 - lambda, lambda / (nc - 1))
      }
      moment <- crossprod(W, W * weights)
      rhs <- crossprod(W, y * weights)
      rhs[1L] <- rhs[1L] + 0.25 * rhs[1L] / moment[1L, 1L]
      drop(E[i, , drop = FALSE] %*%
             solve(moment + diag(0.25, ncol(W)), rhs))
    }, numeric(1L))
    expect_equal(as.vector(H %*% y), oracle, tolerance = 1e-11, info = ukertype)
    expect_equal(as.vector(applied.profile), oracle, tolerance = 1e-11,
                 info = ukertype)
    expect_equal(as.vector(applied.direct), as.vector(H %*% y), tolerance = 1e-11,
                 info = ukertype)
  }
})

test_that("npscoef plot-bootstrap inid helper uses categorical profiles exactly", {
  old <- options(np.messages = FALSE, np.tree = FALSE,
                 np.categorical.compress = FALSE)
  on.exit(options(old), add = TRUE)

  set.seed(20260916L)
  n <- 260L
  B <- 25L
  xdat <- data.frame(x = rnorm(n))
  zdat <- data.frame(
    z = factor(sample(letters[1:4], n, TRUE)),
    o = ordered(sample(1:5, n, TRUE))
  )
  ydat <- 1 + 0.8 * xdat$x + as.numeric(zdat$z) -
    0.4 * as.numeric(zdat$o) + rnorm(n, sd = 0.2)
  bw <- npscoefbw(
    xdat = xdat,
    ydat = ydat,
    zdat = zdat,
    bws = c(0.25, 0.3),
    bandwidth.compute = FALSE,
    regtype = "lc"
  )
  exdat <- xdat[seq_len(31L), , drop = FALSE]
  ezdat <- zdat[seq_len(31L), , drop = FALSE]
  counts <- replicate(B, tabulate(sample.int(n, n, TRUE), nbins = n))
  boot <- getFromNamespace(".np_inid_boot_from_scoef_frozen", "np")

  options(np.tree = FALSE, np.categorical.compress = FALSE)
  dense <- boot(txdat = xdat, ydat = ydat, tzdat = zdat,
                exdat = exdat, ezdat = ezdat, bws = bw, B = B,
                counts = counts, progress.label = "dense")
  options(np.tree = FALSE, np.categorical.compress = TRUE)
  profile <- boot(txdat = xdat, ydat = ydat, tzdat = zdat,
                  exdat = exdat, ezdat = ezdat, bws = bw, B = B,
                  counts = counts, progress.label = "profile")

  expect_true(isTRUE(all.equal(profile$t0, dense$t0, tolerance = 1e-8,
                               check.attributes = FALSE)))
  expect_true(isTRUE(all.equal(profile$t, dense$t, tolerance = 1e-8,
                               check.attributes = FALSE)))
})
