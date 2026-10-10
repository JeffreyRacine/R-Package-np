test_that("LSQ prediction selects retained quantiles and aligned uncertainty", {
  if (!spawn_mpi_slaves(1L)) skip("MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old <- options(np.messages = FALSE); on.exit(options(old), add = TRUE)
  set.seed(101001L)
  d <- data.frame(x = seq(.05, .95, length.out = 36L))
  d$y <- sin(6 * d$x) + rnorm(nrow(d), sd = .15)
  taus <- c(.25, .5, .75)
  model <- nplsqreg(y ~ x, data = d, tau = taus, scale = 1 + d$x,
    delta = .5, bandwidth.compute = FALSE, regtype = "ll", nomad = FALSE)
  # Distinct real scalar fits make an ordering bug observable even when a
  # fixed-delta vector fixture would otherwise have identical columns.
  model$tau.fits <- lapply(seq_along(taus), function(j)
    nplsqreg(bws = .1 + j / 10, txdat = d["x"], tydat = d$y,
      tau = taus[j], scale = 1 + d$x, delta = taus[j],
      bandwidth.compute = FALSE, regtype = "ll", nomad = FALSE))
  nd <- data.frame(x = c(.18, NA, .45, .81), row.names = letters[1:4])
  before <- model
  rng <- .Random.seed
  for (infer in c(FALSE, TRUE)) {
    full <- predict(model, newdata = nd, se.fit = infer)
    expect_identical(predict(model, newdata = nd, se.fit = infer, tau = NULL), full)
    for (ids in list(c(3L, 1L), c(2L, 3L, 1L))) {
      actual <- predict(model, newdata = nd, se.fit = infer, tau = taus[ids])
      if (infer) {
        expect_identical(actual$fit, full$fit[, ids, drop = FALSE])
        expect_identical(actual$se.fit, full$se.fit[, ids, drop = FALSE])
        expect_identical(actual$residual.scale, full$residual.scale[ids])
        expect_identical(actual$df, full$df)
      } else {
        expect_identical(actual, full[, ids, drop = FALSE])
      }
    }
    single <- predict(model, newdata = nd, se.fit = infer, tau = .5)
    reference <- predict(model$tau.fits[[2L]], exdat = nd, se.fit = infer)
    expect_identical(single, reference)
    expect_identical(predict(model, nd, se.fit = infer, tau = .5), single)
    expect_identical(predict(model, newdata = data.frame(wrong = 1),
      exdat = nd, se.fit = infer, tau = .5), single)
    scalar <- model$tau.fits[[2L]]
    expect_identical(predict(scalar, exdat = nd, se.fit = infer, tau = .5),
      predict(scalar, exdat = nd, se.fit = infer))
  }
  expect_identical(model, before)
  expect_identical(.Random.seed, rng)
  expect_identical(predict(unserialize(serialize(model, NULL)), newdata = nd, tau = .75),
    predict(model, newdata = nd, tau = .75))
})

test_that("LSQ prediction rejects invalid or unfitted tau before evaluation", {
  if (!spawn_mpi_slaves(1L)) skip("MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old <- options(np.messages = FALSE); on.exit(options(old), add = TRUE)
  x <- data.frame(x = seq(.05, .95, length.out = 24L))
  model <- nplsqreg(bws = .2, txdat = x, tydat = sin(6 * x$x),
    tau = c(.25, .5, .75), scale = rep(1, nrow(x)), delta = .5,
    bandwidth.compute = FALSE, regtype = "ll", nomad = FALSE)
  for (bad in list(numeric(), "0.5", TRUE, NA_real_, NaN, Inf,
                   0, 1, -.1, 1.1)) {
    expect_error(predict(model, tau = bad,
      newdata = stop("newdata was forced")), "finite numeric values")
  }
  expect_error(predict(model, tau = c(.5, .5),
    newdata = stop("newdata was forced")), "duplicate 'tau'")
  for (bad in list(.4, c(.25, .4), .5 + 1e-12)) {
    expect_error(predict(model, tau = bad,
      newdata = stop("newdata was forced")), "was not fitted")
  }
  expect_error(predict(model$tau.fits[[1L]], tau = .5), "was not fitted")
  expect_error(predict(model, tau = .5, total_nonsense = TRUE), "unused")
  expect_error(predict(model, tau = .5, se.fit = FALSE, se = TRUE), "uses se.fit")
  expect_error(predict(model, tau = .5, se.fit = NA), "se.fit")
})

test_that("selected LSQ children retain parent formula and omission ownership", {
  if (!spawn_mpi_slaves(1L)) skip("MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  old <- options(np.messages = FALSE); on.exit(options(old), add = TRUE)
  set.seed(101002L)
  d <- data.frame(x = seq(.05, .95, length.out = 36L), y = rnorm(36L))
  d$x[3L] <- NA_real_
  counter <- new.env(parent = emptyenv()); counter$n <- 0L
  counted <- function(x) { counter$n <- counter$n + 1L; x }
  model <- nplsqreg(y ~ counted(x), data = d, tau = c(.25, .5, .75),
    scale = rep(1, nrow(d)), delta = .5, bandwidth.compute = FALSE,
    regtype = "ll", nomad = FALSE, na.action = na.exclude)
  nd <- data.frame(x = c(.2, NA, .5, .8), row.names = LETTERS[1:4])
  for (infer in c(FALSE, TRUE)) {
    for (request in list(.5, c(.75, .25))) {
      counter$n <- 0L
      actual <- predict(model, newdata = nd, se.fit = infer, tau = request)
      expect_identical(counter$n, 1L)
      value <- if (infer) actual$fit else actual
      expect_identical(NROW(value), 4L)
      expect_true(all(is.na(if (is.matrix(value)) value[2L, ] else value[2L])))
      training <- predict(model, se.fit = infer, tau = request)
      value <- if (infer) training$fit else training
      expect_identical(NROW(value), nrow(d))
      expect_true(all(is.na(if (is.matrix(value)) value[3L, ] else value[3L])))
    }
  }
  counter$n <- 0L
  expect_error(predict(model, newdata = nd, tau = .4), "was not fitted")
  expect_identical(counter$n, 0L)
  broken <- model; broken$tau.fits <- NULL
  expect_error(predict(broken, newdata = nd, tau = .5), "per-tau fit state")
})
