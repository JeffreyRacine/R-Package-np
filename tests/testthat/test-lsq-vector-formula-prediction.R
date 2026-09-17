test_that("vector least-squares quantile prediction prepares newdata once", {
  if (!spawn_mpi_slaves(2L)) skip("MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old <- options(np.messages = FALSE); on.exit(options(old), add = TRUE)
  set.seed(1010L)
  d <- data.frame(x = runif(40), y = rnorm(40))
  counts <- new.env(parent = emptyenv()); counts$n <- 0L
  counted <- function(x) { counts$n <- counts$n + 1L; x }
  f <- y ~ counted(x)
  model <- nplsqreg(f, data = d, tau = c(.25, .75), scale = rep(1, nrow(d)),
    delta = .5, bandwidth.compute = FALSE, regtype = "ll", nomad = FALSE)
  expect_identical(counts$n, 1L)
  nd <- data.frame(x = seq(.15, .85, length.out = 9))
  ex <- nd; names(ex) <- "counted(x)"
  for (infer in c(FALSE, TRUE)) {
    counts$n <- 0L; rng <- .Random.seed
    actual <- predict(model, newdata = nd, se.fit = infer)
    expect_identical(counts$n, 1L)
    expect_identical(.Random.seed, rng)
    references <- lapply(model$tau.fits, predict, exdat = ex, se.fit = infer)
    expected <- if (infer) do.call(cbind, lapply(references, `[[`, "fit")) else do.call(cbind, references)
    colnames(expected) <- colnames(if (infer) actual$fit else actual)
    expect_identical(if (infer) actual$fit else actual, expected)
    if (infer) {
      expected.se <- do.call(cbind, lapply(references, `[[`, "se.fit"))
      colnames(expected.se) <- colnames(actual$se.fit)
      expect_identical(actual$se.fit, expected.se)
    }
    counts$n <- 0L
    explicit <- predict(model, newdata = data.frame(wrong = 1:9), exdat = ex, se.fit = infer)
    expect_identical(counts$n, 0L)
    expect_identical(explicit, actual)
  }
  counts$n <- 0L
  expect_error(predict(model, newdata = nd, se.fit = FALSE, se = TRUE), "uses se.fit")
  expect_identical(counts$n, 0L)
  expect_error(predict(model, newdata = data.frame(wrong = 1:9)), "columns")
  broken <- model; broken$tau.fits <- NULL
  expect_error(predict(broken, newdata = nd), "per-tau fit state")
})

test_that("vector least-squares prediction shares stochastic evaluation values", {
  if (!spawn_mpi_slaves(2L)) skip("MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old <- options(np.messages = FALSE); on.exit(options(old), add = TRUE)
  set.seed(1011L)
  d <- data.frame(x = runif(40), y = rnorm(40))
  jittered <- function(x) x + rnorm(length(x), sd = .01)
  model <- nplsqreg(y ~ jittered(x), data = d, tau = c(.25, .75),
    scale = rep(1, nrow(d)), delta = .5, bandwidth.compute = FALSE,
    regtype = "ll", nomad = FALSE)
  nd <- data.frame(x = seq(.15, .85, length.out = 9))
  set.seed(1012L)
  ex <- data.frame(jittered(nd$x)); names(ex) <- "jittered(x)"
  expected.rng <- .Random.seed
  reference <- predict(model, exdat = ex, se.fit = TRUE)
  set.seed(1012L)
  actual <- predict(model, newdata = nd, se.fit = TRUE)
  expect_identical(.Random.seed, expected.rng)
  expect_identical(actual, reference)
})

test_that("vector prediction retains mixed time-indexed and native input semantics", {
  skip_on_cran()
  spawn_mpi_slaves(2L)
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old <- options(np.messages = FALSE); on.exit(options(old), add = TRUE)
  set.seed(1013L)
  d <- data.frame(x = runif(48), v = runif(48),
                  u = factor(rep(c("a", "b"), 24)), y = rnorm(48))
  model <- nplsqreg(y ~ x + v + u, data = d, tau = c(.25, .75),
    scale = rep(1, nrow(d)), delta = .5, bandwidth.compute = FALSE,
    regtype = "ll", nomad = FALSE)
  ts.model <- nplsqreg(y ~ x + v, data = d, tau = c(.25, .75),
    scale = rep(1, nrow(d)), delta = .5, bandwidth.compute = FALSE,
    regtype = "ll", nomad = FALSE)
  nd <- data.frame(x = ts(seq(.1, .9, length.out = 12), start = 1),
                   v = ts(seq(.2, .8, length.out = 12), start = 3))
  ex <- data.frame(x = as.numeric(nd$x)[3:12], v = as.numeric(nd$v)[1:10],
                   row.names = NULL)
  expect_identical(predict(ts.model, newdata = nd, se.fit = TRUE),
                   predict(ts.model, exdat = ex, se.fit = TRUE))
  ex <- d[1:10, c("x", "v", "u")]
  native <- nplsqreg(bws = model$bws, txdat = d[c("x", "v", "u")],
    tydat = d$y)
  native$bws$formula <- NULL
  for (j in seq_along(native$tau.fits))
    native$tau.fits[[j]]$bws$formula <- NULL
  expect_identical(predict(native, newdata = ex), predict(native, exdat = ex))
  ex$x[2] <- NA_real_
  expected <- lapply(model$tau.fits, predict, exdat = ex)
  actual <- predict(model, newdata = ex)
  expected <- do.call(cbind, expected); colnames(expected) <- colnames(actual)
  expect_identical(actual, expected)
  expect_identical(dim(actual), c(9L, 2L))
})
