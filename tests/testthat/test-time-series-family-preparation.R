test_that("kernel sums and copulas use the aligned joint sample", {
  skip_if_not(spawn_mpi_slaves())
  on.exit(close_mpi_slaves(), add = TRUE)
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  set.seed(42)
  y <- ts(rnorm(28), frequency = 4)
  x <- data.frame(a = as.numeric(y)[2:27], b = as.numeric(y)[1:26])
  response <- as.numeric(y)[3:28]
  h <- c(.8, .8)
  direct <- npksum(y ~ lag(y, -1) + lag(y, -2), bws = h)
  control <- npksum(txdat = x, tydat = response, bws = h)
  expect_identical(direct$ksum, control$ksum)
  expect_length(direct$ksum, 26L)
  keep <- seq_len(26) > 2L
  direct <- npksum(y ~ lag(y, -1) + lag(y, -2), bws = h, subset = keep)
  control <- npksum(txdat = x[keep, ], tydat = response[keep], bws = h)
  expect_identical(direct$ksum, control$ksum)
  joint <- data.frame(y = response, b = x$b)
  dat <- data.frame(y = as.numeric(y))
  dat$y <- y
  for (target in c("distribution", "density")) {
    bw.fun <- if (target == "distribution") npudistbw else npudensbw
    bw <- bw.fun(dat = joint, bws = h, bandwidth.compute = FALSE)
    control <- npcopula(bws = bw, data = joint, target = target, evaluation = "sample", se = TRUE)
    direct <- npcopula(~y + lag(y, -2), data = dat, bws = h, bandwidth.compute = FALSE,
      target = target, evaluation = "sample", se = TRUE)
    expect_identical(direct$ntrain, 26L)
    expect_equal(direct$copula, control$copula, tolerance = 1e-14)
    expect_length(direct$copulaerr, 26L)
    expect_equal(direct$copulaerr, control$copulaerr, tolerance = 1e-14)
  }
})

test_that("least-squares quantile formulas retain aligned training and evaluation", {
  skip_if_not(spawn_mpi_slaves())
  on.exit(close_mpi_slaves(), add = TRUE)
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  set.seed(42)
  y <- ts(rnorm(28), frequency = 4)
  f <- y ~ lag(y, -1) + lag(y, -2)
  x <- data.frame(a = as.numeric(y)[2:27], b = as.numeric(y)[1:26])
  response <- as.numeric(y)[3:28]
  names(x) <- c("lag(y, -1)", "lag(y, -2)")
  bw <- nplsqregbw(f, nmulti = 1L, scale = rep(1, 26), optim.control = list(maxit = 2L))
  expect_identical(unname(as.matrix(bw$xdat)), unname(as.matrix(x)))
  expect_identical(as.numeric(bw$ydat), response)
  fit <- nplsqreg(f, nmulti = 1L, scale = rep(1, 26), optim.control = list(maxit = 2L))
  control <- nplsqreg(bws = fit$bws, txdat = x, tydat = response)
  expect_equal(fitted(fit), fitted(control), tolerance = 1e-14)
  new <- data.frame(y = as.numeric(y)[1:18])
  new$y <- ts(new$y, frequency = 4, start = c(2001, 1))
  ex <- data.frame(a = as.numeric(new$y)[2:18], b = as.numeric(new$y)[1:17])
  names(ex) <- names(x)
  control <- nplsqreg(bws = fit$bws, txdat = x, tydat = response, exdat = ex)
  expect_equal(predict(fit, newdata = new), fitted(control), tolerance = 1e-14)
})

test_that("both IV formula preparers intersect before separating roles", {
  parse <- getFromNamespace(".np_iv_parse_formula", "npRmpi")
  frame <- getFromNamespace(".np_iv_training_frame", "npRmpi")
  eval.frame <- getFromNamespace(".np_iv_eval_frame", "npRmpi")
  set.seed(42)
  y <- ts(rnorm(28), frequency = 4)
  for (owner in c("npregiv()", "npregivderiv()")) {
    parsed <- parse(y ~ lag(y, -1) | lag(y, -2), owner)
    actual <- frame(parsed, quote(f()), NULL, owner, environment())
    expect_identical(unname(as.matrix(actual)), cbind(as.numeric(y)[3:28],
      as.numeric(y)[2:27], as.numeric(y)[1:26]))
    new <- data.frame(y = rnorm(18))
    new$y <- ts(new$y, frequency = 4, start = c(2001, 1))
    actual <- eval.frame(parsed, new, names(parsed$roles), NULL, owner)
    expect_identical(unname(as.matrix(actual)), cbind(as.numeric(new$y)[2:18],
      as.numeric(new$y)[1:17]))
  }
})
