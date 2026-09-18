test_that("smooth-coefficient formula transactions prepare expressions once", {
  old <- options(np.messages = FALSE); on.exit(options(old), add = TRUE)
  set.seed(974L)
  d <- data.frame(x = runif(40), z = runif(40), y = rnorm(40))
  counts <- new.env(parent = emptyenv()); counts$n <- 0L
  counted <- function(x) { counts$n <- counts$n + 1L; x }
  for (three in c(FALSE, TRUE)) {
    f <- if (three) y ~ counted(x) | z else y ~ counted(x)
    counts$n <- 0L
    bw <- npscoefbw(f, data = d, bws = .5, bandwidth.compute = FALSE)
    expect_identical(counts$n, 1L)
    expect_false(grepl("getFromNamespace", paste(deparse(attr(bw$terms, "predvars")), collapse = "")))
    counts$n <- 0L
    ref <- npscoef(bw, se = TRUE, betas = TRUE)
    expect_identical(counts$n, 0L)
    for (args in list(list(f, data = d, bws = .5),
        list(formula = f, data = d, bws = .5),
        list(data = d, bws = .5, formula = f))) {
      counts$n <- 0L; rng <- .Random.seed
      actual <- do.call(npscoef, c(args, list(bandwidth.compute = FALSE, se = TRUE, betas = TRUE)))
      expect_identical(counts$n, 1L)
      expect_identical(.Random.seed, rng)
      for (field in c("mean", "merr", "beta", "grad", "gerr", "R2", "MSE"))
        expect_identical(actual[[field]], ref[[field]], info = field)
      expect_false(grepl("formula.state", paste(deparse(actual$bws$call), collapse = "")))
    }
    counts$n <- 0L
    evaluated <- npscoef(bw, newdata = d[1:7, ], se = TRUE)
    expect_identical(counts$n, 1L)
    counts$n <- 0L
    prediction <- predict(ref, newdata = d[1:7, ], se.fit = TRUE)
    expect_identical(counts$n, 1L)
    expect_identical(prediction$fit, fitted(evaluated))
    expect_identical(prediction$se.fit, se(evaluated))
    expect_error(predict(ref, newdata = data.frame(wrong = 1:7)), "columns.*x")
  }
})

test_that("smooth-coefficient formulas retain the indexed sample and trained terms", {
  old <- options(np.messages = FALSE); on.exit(options(old), add = TRUE)
  set.seed(975L)
  y <- ts(rnorm(42), frequency = 4)
  f <- y ~ lag(y, -1) | lag(y, -2)
  x <- data.frame(as.numeric(y)[2:41]); names(x) <- "lag(y, -1)"
  z <- data.frame(as.numeric(y)[1:40]); names(z) <- "lag(y, -2)"
  response <- as.numeric(y)[3:42]
  nd <- data.frame(y = rnorm(18)); nd$y <- ts(nd$y, frequency = 4)
  ex <- data.frame(as.numeric(nd$y)[2:18]); names(ex) <- names(x)
  ez <- data.frame(as.numeric(nd$y)[1:17]); names(ez) <- names(z)
  for (type in c("fixed", "generalized_nn", "adaptive_nn")) {
    h <- if (type == "fixed") .8 else 20
    bw <- npscoefbw(f, bws = h, bandwidth.compute = FALSE, bwtype = type)
    fit <- npscoef(bw, se = TRUE, betas = TRUE)
    ref <- npscoef(bw, txdat = x, tydat = response, tzdat = z, se = TRUE, betas = TRUE)
    value <- npscoef(bw, newdata = nd, se = TRUE, betas = TRUE)
    oracle <- npscoef(bw, txdat = x, tydat = response, tzdat = z,
      exdat = ex, ezdat = ez, se = TRUE, betas = TRUE)
    for (field in c("mean", "merr", "beta", "grad", "gerr", "R2", "MSE")) {
      expect_equal(fit[[field]], ref[[field]], tolerance = 1e-12, info = paste(type, field))
      expect_equal(value[[field]], oracle[[field]], tolerance = 1e-12, info = paste(type, field))
    }
    expect_length(fitted(fit), 40L)
    expect_length(fitted(value), 17L)
  }
  d <- data.frame(y = rnorm(36), x = runif(36), z = runif(36))
  f <- y ~ poly(x, degree = 1) | z
  bw <- npscoefbw(f, data = d, bws = .5, bandwidth.compute = FALSE)
  mf <- model.frame(y ~ poly(x, degree = 1) + z, data = d)
  # Retained training uses the initial poly QR values; only newdata uses the
  # trained recurrence. The native oracle must describe those same two samples.
  new <- transform(d[1:8, ], x = x + .1)
  em <- model.frame(delete.response(attr(mf, "terms")), data = new)
  value <- npscoef(bw, newdata = new, se = FALSE)
  oracle <- npscoef(bw, txdat = mf[2], tydat = mf[[1]], tzdat = mf[3],
    exdat = em[1], ezdat = em[2], se = FALSE)
  expect_identical(fitted(value), fitted(oracle))
})

test_that("smooth-coefficient one-call formulas use the first stochastic sample", {
  old <- options(np.messages = FALSE); on.exit(options(old), add = TRUE)
  set.seed(976L)
  d <- data.frame(x = runif(36), z = runif(36), y = rnorm(36))
  jittered <- function(x) x + rnorm(length(x), sd = .01)
  f <- y ~ jittered(x) | z
  set.seed(977L)
  x <- data.frame(jittered(d$x)); names(x) <- "jittered(x)"
  rng <- .Random.seed
  ref <- npscoef(txdat = x, tydat = d$y, tzdat = d["z"],
    bws = .5, bandwidth.compute = FALSE, se = FALSE)
  set.seed(977L)
  fit <- npscoef(f, data = d, bws = .5, bandwidth.compute = FALSE, se = FALSE)
  expect_identical(.Random.seed, rng)
  expect_identical(fitted(fit), fitted(ref))
})
