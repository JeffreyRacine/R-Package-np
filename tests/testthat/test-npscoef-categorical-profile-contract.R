npscoef_profile_oracle <- function(xdat, ydat, zdat, bws, exdat = NULL,
                                   ezdat = NULL) {
  NZD <- getFromNamespace("NZD", "npRmpi")
  npRidgeSequenceAdditive <- getFromNamespace("npRidgeSequenceAdditive", "npRmpi")
  .np_cat_profile_code_matrix <-
    getFromNamespace(".np_cat_profile_code_matrix", "npRmpi")
  .np_cat_profile_keys <-
    getFromNamespace(".np_cat_profile_keys", "npRmpi")
  .np_regression_cat_profile_kernel_matrix <-
    getFromNamespace(".np_regression_cat_profile_kernel_matrix", "npRmpi")
  .np_cat_profile_rowsum <-
    getFromNamespace(".np_cat_profile_rowsum", "npRmpi")

  if (is.null(exdat))
    exdat <- xdat
  if (is.null(ezdat))
    ezdat <- zdat

  W.train <- as.matrix(data.frame(1, xdat))
  W.eval <- as.matrix(data.frame(1, exdat))
  yW <- cbind(ydat, W.train)
  train.codes <- .np_cat_profile_code_matrix(zdat)
  eval.codes <- .np_cat_profile_code_matrix(ezdat)
  train.keys <- .np_cat_profile_keys(train.codes)
  profile.keys <- unique(train.keys)
  train.id <- match(train.keys, profile.keys)
  train.rep <- match(profile.keys, train.keys)
  G <- length(profile.keys)
  L <- .np_regression_cat_profile_kernel_matrix(
    eval.codes = eval.codes,
    train.codes = train.codes[train.rep, , drop = FALSE],
    xdat = zdat[train.rep, , drop = FALSE],
    bws = bws
  )

  p <- ncol(yW)
  cross.profile <- matrix(0.0, nrow = G, ncol = p * p)
  for (j in seq_len(p)) {
    for (k in seq_len(p)) {
      cross.profile[, (j - 1L) * p + k] <-
        .np_cat_profile_rowsum(yW[, j] * yW[, k], train.id, G)[, 1L]
    }
  }
  main.flat <- L %*% cross.profile
  main <- array(t(main.flat), dim = c(p, p, nrow(ezdat)))
  tyw <- main[-1L, 1L, , drop = FALSE]
  dim(tyw) <- c(dim(tyw)[1L], dim(tyw)[3L])
  tww <- main[-1L, -1L, , drop = FALSE]

  ridge.grid <- npRidgeSequenceAdditive(n.train = nrow(xdat), cap = 1.0)
  coef <- matrix(NA_real_, nrow = ncol(W.eval), ncol = nrow(W.eval))
  for (i in seq_len(nrow(W.eval))) {
    for (ridge in ridge.grid) {
      ridge.val <- ridge * tyw[1L, i] / NZD(tww[1L, 1L, i])
      theta <- tryCatch(
        solve(tww[, , i] + diag(rep(ridge, nrow(tyw))),
              tyw[, i] + c(ridge.val, rep(0, nrow(tyw) - 1L))),
        error = function(e) e
      )
      if (!inherits(theta, "error")) {
        coef[, i] <- theta
        break
      }
    }
  }
  list(mean = as.vector(colSums(t(W.eval) * coef)),
       beta = t(coef))
}

test_that("npscoef all-categorical profile route preserves fitted values", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  options(npRmpi.autodispatch = FALSE)
  on.exit(close_mpi_slaves(), add = TRUE)
  old <- options(np.messages = FALSE, np.tree = FALSE,
                 np.categorical.compress = TRUE)
  on.exit(options(old), add = TRUE)

  set.seed(20260516L)
  n <- 100L
  xdat <- data.frame(x = rnorm(n))
  zdat <- data.frame(
    z = factor(sample(letters[1:3], n, TRUE)),
    o = ordered(sample(1:4, n, TRUE))
  )
  ydat <- 1 + 0.7 * xdat$x +
    rowSums(as.data.frame(lapply(zdat, as.numeric))) +
    rnorm(n, sd = 0.15)
  bw <- npscoefbw(
    xdat = xdat,
    ydat = ydat,
    zdat = zdat,
    bws = c(0.25, 0.3),
    bandwidth.compute = FALSE,
    regtype = "lc"
  )
  fit <- npscoef(bws = bw, txdat = xdat, tydat = ydat, tzdat = zdat,
                 errors = FALSE, iterate = FALSE, betas = TRUE)
  oracle <- npscoef_profile_oracle(xdat, ydat, zdat, bw)
  expect_equal(fit$mean, oracle$mean, tolerance = 1e-8)
  expect_equal(fit$beta, oracle$beta, tolerance = 1e-8)
})

test_that("npscoef all-categorical profile route preserves evaluation values", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  options(npRmpi.autodispatch = FALSE)
  on.exit(close_mpi_slaves(), add = TRUE)
  old <- options(np.messages = FALSE, np.tree = FALSE,
                 np.categorical.compress = TRUE)
  on.exit(options(old), add = TRUE)

  set.seed(20260616L)
  n <- 100L
  xdat <- data.frame(x = rnorm(n))
  zdat <- data.frame(o = ordered(sample(1:5, n, TRUE)))
  ydat <- 1 + 0.7 * xdat$x + as.numeric(zdat$o) + rnorm(n, sd = 0.15)
  bw <- npscoefbw(
    xdat = xdat,
    ydat = ydat,
    zdat = zdat,
    bws = 0.3,
    bandwidth.compute = FALSE,
    regtype = "lc",
    okertype = "racineliyan"
  )
  exdat <- data.frame(x = seq(min(xdat$x), max(xdat$x), length.out = 31L))
  ezdat <- zdat[seq_len(31L), , drop = FALSE]
  fit <- npscoef(bws = bw, txdat = xdat, tydat = ydat, tzdat = zdat,
                 exdat = exdat, ezdat = ezdat,
                 errors = FALSE, iterate = FALSE, betas = TRUE)
  oracle <- npscoef_profile_oracle(xdat, ydat, zdat, bw,
                                   exdat = exdat, ezdat = ezdat)
  expect_equal(fit$mean, oracle$mean, tolerance = 1e-8)
  expect_equal(fit$beta, oracle$beta, tolerance = 1e-8)
})

test_that("npscoef all-categorical profile route preserves bandwidth CV", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  options(npRmpi.autodispatch = TRUE)
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old <- options(np.messages = FALSE, np.tree = FALSE,
                 np.categorical.compress = FALSE)
  on.exit(options(old), add = TRUE)

  set.seed(20261017L)
  n <- 80L
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

test_that("npscoefhat apply uses categorical profile route", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  options(npRmpi.autodispatch = FALSE)
  on.exit(close_mpi_slaves(), add = TRUE)
  old <- options(np.messages = FALSE, np.tree = FALSE,
                 np.categorical.compress = TRUE)
  on.exit(options(old), add = TRUE)

  set.seed(20260716L)
  n <- 100L
  xdat <- data.frame(x = rnorm(n))
  zdat <- data.frame(
    z = factor(sample(letters[1:3], n, TRUE)),
    o = ordered(sample(1:4, n, TRUE))
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
  got <- npscoefhat(bws = bw, txdat = xdat, tzdat = zdat,
                   y = rhs, output = "apply")
  oracle <- npscoef_profile_oracle(xdat, ydat, zdat, bw)
  expect_equal(got[, 1L], oracle$mean, tolerance = 1e-8)
})

test_that("npscoef plot-bootstrap inid helper uses categorical profiles exactly", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  options(npRmpi.autodispatch = FALSE)
  on.exit(close_mpi_slaves(), add = TRUE)
  old <- options(np.messages = FALSE, np.tree = FALSE,
                 np.categorical.compress = FALSE)
  on.exit(options(old), add = TRUE)

  set.seed(20260816L)
  n <- 80L
  B <- 8L
  xdat <- data.frame(x = rnorm(n))
  zdat <- data.frame(
    z = factor(sample(letters[1:3], n, TRUE)),
    o = ordered(sample(1:4, n, TRUE))
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
  exdat <- xdat[seq_len(19L), , drop = FALSE]
  ezdat <- zdat[seq_len(19L), , drop = FALSE]
  counts <- replicate(B, tabulate(sample.int(n, n, TRUE), nbins = n))
  boot <- getFromNamespace(".np_inid_boot_from_scoef_frozen", "npRmpi")

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
