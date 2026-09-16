test_that("smooth-coefficient automatic formulas search on one common sample", {
  old <- options(np.messages = FALSE); on.exit(options(old), add = TRUE)
  set.seed(978L)
  d <- data.frame(x = runif(24), z = runif(24), y = rnorm(24))
  counts <- new.env(parent = emptyenv()); counts$n <- 0L
  counted <- function(x) { counts$n <- counts$n + 1L; x }
  f <- y ~ counted(x) | z
  results <- list()
  for (named in c(FALSE, TRUE)) {
    counts$n <- 0L; set.seed(979L)
    args <- if (named) list(data = d, formula = f) else list(f, data = d)
    value <- do.call(npscoef, c(args, list(nmulti = 1L, optim.maxit = 5L, se = FALSE)))
    expect_identical(counts$n, 1L)
    results[[length(results) + 1L]] <- list(result = value$mean,
      bw = value$bws$bw, fval = value$bws$fval,
      feval = value$bws$num.feval, rng = .Random.seed)
  }
  expect_identical(results[[1L]], results[[2L]])
})

test_that("smooth-coefficient subsets and omission metadata retain their sample", {
  old <- options(np.messages = FALSE); on.exit(options(old), add = TRUE)
  set.seed(980L)
  d <- data.frame(x = runif(40), z = runif(40), y = rnorm(40))
  d$y[7L] <- NA_real_
  bw <- npscoefbw(y ~ x | z, data = d, subset = x > .2, na.action = na.omit,
    bws = .5, bandwidth.compute = FALSE)
  mf <- model.frame(y ~ x + z, data = d, subset = x > .2, na.action = na.omit)
  fit <- npscoef(bw, se = TRUE)
  ref <- npscoef(bw, txdat = mf["x"], tydat = mf$y, tzdat = mf["z"], se = TRUE)
  expect_identical(fitted(fit), fitted(ref))
  expect_identical(se(fit), se(ref))
  expect_identical(bw$rows.omit, as.vector(attr(mf, "na.action")))
  expect_identical(fit$ntrain, nrow(mf))
  expect_error(npscoef(formula = y ~ missing.symbol | z, data = d,
    bws = .5, bandwidth.compute = FALSE), "missing.symbol")
  good <- npscoef(bw, txdat = mf["x"], tydat = mf$y, tzdat = mf["z"], se = FALSE)
  expect_identical(fitted(good), fitted(ref))
})
