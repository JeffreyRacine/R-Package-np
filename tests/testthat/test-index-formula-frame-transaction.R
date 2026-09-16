test_that("single-index formulas prepare one sample independently of argument order", {
  old <- options(np.messages = FALSE); on.exit(options(old), add = TRUE)
  set.seed(925L)
  d <- data.frame(x = runif(30), z = rnorm(30), y = rnorm(30))
  counter <- new.env(parent = emptyenv()); counter$n <- 0L
  counted <- function(x) { counter$n <- counter$n + 1L; x }
  f <- y ~ counted(x) + z
  h <- c(1, .2, .8)
  bw <- npindexbw(f, data = d, bws = h, bandwidth.compute = FALSE)
  expect_identical(counter$n, 1L)
  expect_false(grepl("getFromNamespace", paste(deparse(attr(bw$terms, "predvars")), collapse = "")))
  counter$n <- 0L
  fit <- npindex(bw, se = FALSE)
  expect_identical(counter$n, 1L)
  counter$n <- 0L
  evaluated <- npindex(bw, newdata = d[1:5, ], se = FALSE)
  expect_identical(counter$n, 2L)
  expect_equal(predict(fit, newdata = d[1:5, ]), fitted(evaluated), tolerance = 0)
  for (args in list(list(f, data = d, bws = h),
      list(formula = f, data = d, bws = h),
      list(data = d, formula = f, bws = h))) {
    counter$n <- 0L; rng <- .Random.seed
    value <- do.call(npindex, c(args, list(se = FALSE)))
    expect_identical(counter$n, 1L)
    expect_identical(.Random.seed, rng)
    expect_identical(value$beta, fit$beta)
    expect_identical(value$bw, fit$bw)
    expect_identical(fitted(value), fitted(fit))
  }
  output <- list()
  for (named in c(FALSE, TRUE)) {
    counter$n <- 0L; set.seed(926L)
    args <- if (named) list(data = d, formula = f) else list(f, data = d)
    output[[length(output) + 1L]] <- do.call(npindex,
      c(args, list(nmulti = 1L, itmax = 8L, se = FALSE)))
    expect_identical(counter$n, 1L)
    output[[length(output)]]$rng <- .Random.seed
  }
  for (field in c("beta", "bw", "mean", "rng"))
    expect_identical(output[[1L]][[field]], output[[2L]][[field]], info = field)
  expect_identical(output[[1L]]$bws$fval, output[[2L]]$bws$fval)
  expect_identical(output[[1L]]$bws$num.feval, output[[2L]]$bws$num.feval)
})

test_that("single-index training and evaluation agree with integer time indices", {
  old <- options(np.messages = FALSE); on.exit(options(old), add = TRUE)
  set.seed(927L)
  y <- ts(rnorm(36), frequency = 4)
  f <- y ~ lag(y, -1) + lag(y, -2)
  x <- data.frame(a = as.numeric(y)[2:35], b = as.numeric(y)[1:34])
  names(x) <- attr(terms(f), "term.labels")
  response <- as.numeric(y)[3:36]
  new <- data.frame(y = rnorm(20)); new$y <- ts(new$y, frequency = 4)
  ex <- data.frame(a = as.numeric(new$y)[2:20], b = as.numeric(new$y)[1:19])
  names(ex) <- names(x)
  for (type in c("fixed", "generalized_nn", "adaptive_nn")) {
    h <- c(1, .2, if (type == "fixed") .8 else 16)
    bw <- npindexbw(f, bws = h, bandwidth.compute = FALSE, bwtype = type)
    fit <- npindex(bw, se = TRUE, gradients = TRUE)
    native <- npindex(bw, txdat = x, tydat = response, se = TRUE, gradients = TRUE)
    evaluated <- npindex(bw, newdata = new, se = TRUE, gradients = TRUE)
    oracle <- npindex(bw, txdat = x, tydat = response, exdat = ex,
      se = TRUE, gradients = TRUE)
    for (field in c("mean", "merr", "grad", "gerr", "betavcov", "R2", "MSE")) {
      expect_equal(fit[[field]], native[[field]], tolerance = 1e-12, info = paste(type, field))
      expect_equal(evaluated[[field]], oracle[[field]], tolerance = 1e-12, info = paste(type, field))
    }
    expect_length(fitted(fit), 34L)
    expect_length(fitted(evaluated), 19L)
    expect_equal(predict(fit, newdata = new), fitted(oracle), tolerance = 1e-12)
    expect_error(predict(fit, newdata = data.frame(other = 1:20)), "columns.*y")
  }
})

test_that("single-index retained transform metadata is used for prediction", {
  old <- options(np.messages = FALSE); on.exit(options(old), add = TRUE)
  set.seed(928L)
  d <- data.frame(x = runif(32), z = rnorm(32), y = rnorm(32))
  f <- y ~ poly(x, degree = 1) + z
  h <- c(1, .2, .8)
  bw <- npindexbw(f, data = d, bws = h, bandwidth.compute = FALSE)
  mf <- model.frame(f, data = d)
  nd <- d[1:6, ]; nd$x <- nd$x + .1
  em <- model.frame(delete.response(attr(mf, "terms")), data = nd)
  fit <- npindex(bw, newdata = nd, se = FALSE)
  native <- npindex(bw, txdat = mf[-1], tydat = mf[[1]], exdat = em, se = FALSE)
  expect_identical(fitted(fit), fitted(native))
  direct <- npindex(f, data = d, bws = h, se = FALSE)
  expect_identical(predict(direct, newdata = nd), fitted(native))
  dot <- npindex(y ~ ., data = d, bws = h, se = FALSE)
  ordinary <- npindex(y ~ x + z, data = d, bws = h, se = FALSE)
  expect_identical(fitted(dot), fitted(ordinary))
})
