test_that("npksum prepares formula environments and mixed series before transport", {
  skip_if_not(spawn_mpi_slaves())
  on.exit(close_mpi_slaves(), add = TRUE)
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  set.seed(821L)
  y <- ts(rnorm(30), start = c(2000, 1), frequency = 4)
  x <- ts(rnorm(30), start = c(2000, 1), frequency = 4)
  f <- factor(rep(c("a", "b"), length.out = 29))
  oracle.x <- data.frame(x = as.numeric(x)[1:29])
  oracle.y <- as.numeric(y)[2:30]
  reference <- npksum(txdat = oracle.x, tydat = oracle.y, bws = .8)
  from.env <- npksum(y ~ lag(x, -1), data = environment(), bws = .8)
  expect_equal(from.env$ksum, reference$ksum, tolerance = 1e-14)
  from.list <- npksum(y ~ lag(x, -1), data = list(y = y, x = x), bws = .8)
  expect_equal(from.list$ksum, reference$ksum, tolerance = 1e-14)
  from.local <- local({
    response <- y
    predictor <- x
    npksum(response ~ lag(predictor, -1), bws = .8)
  })
  expect_equal(from.local$ksum, reference$ksum, tolerance = 1e-14)
  mixed <- npksum(y ~ lag(x, -1) + f, bws = c(.8, .3))
  mixed.reference <- npksum(txdat = data.frame(x = oracle.x$x, f = f),
    tydat = oracle.y, bws = c(.8, .3))
  expect_equal(mixed$ksum, mixed.reference$ksum, tolerance = 1e-14)
  expect_length(mixed$ksum, 29L)
  dat <- data.frame(y = as.numeric(y), x = as.numeric(x))
  dat$x <- x
  dat$y <- y
  from.frame <- npksum(y ~ lag(x, -1), data = dat, bws = .8)
  expect_equal(from.frame$ksum, reference$ksum, tolerance = 1e-14)
  keep <- seq_len(29) > 3L
  subset <- npksum(y ~ lag(x, -1), data = environment(), subset = keep, bws = .8)
  expect_equal(subset$ksum,
    npksum(txdat = oracle.x[keep, , drop = FALSE], tydat = oracle.y[keep], bws = .8)$ksum,
    tolerance = 1e-14)
})

test_that("npksum formula failures leave the pool usable", {
  skip_if_not(spawn_mpi_slaves())
  on.exit(close_mpi_slaves(), add = TRUE)
  y <- ts(seq_len(12), frequency = 4)
  x <- ts(cbind(a = 101:112, b = 201:212), frequency = 4)
  expect_error(npksum(y ~ .np.missing.a4.symbol, bws = .8),
               "object '.np.missing.a4.symbol' not found", fixed = TRUE)
  expect_warning(expect_error(npksum(y ~ lag(x, -1), bws = .8),
    "do not match|continuous variable|invalid"), NA)
  healthy <- npksum(y ~ lag(y, -1), bws = .8)
  reference <- npksum(txdat = data.frame(x = as.numeric(y)[1:11]),
    tydat = as.numeric(y)[2:12], bws = .8)
  expect_equal(healthy$ksum, reference$ksum, tolerance = 1e-14)
})

test_that("npksum formula leaves retain operators and bandwidth topology", {
  skip_if_not(spawn_mpi_slaves())
  on.exit(close_mpi_slaves(), add = TRUE)
  set.seed(830L)
  dat <- data.frame(x = runif(40), y = rnorm(40))
  ex <- data.frame(x = c(.18, .47, .81))
  for (type in c("fixed", "generalized_nn", "adaptive_nn")) {
    bw <- if (type == "fixed") .2 else 8
    for (operator in c("normal", "derivative", "integral")) {
      direct <- npksum(y ~ x, data = dat, newdata = ex,
        bws = bw, bwtype = type, operator = operator, weights = matrix(seq_len(40) / 40))
      reference <- npksum(txdat = dat["x"], tydat = dat$y, exdat = ex,
        bws = bw, bwtype = type, operator = operator, weights = matrix(seq_len(40) / 40))
      expect_equal(direct$ksum, reference$ksum, tolerance = 1e-14)
    }
    direct <- npksum(y ~ log(x + 1), data = dat, bws = bw,
      bwtype = type, leave.one.out = TRUE)
    reference <- npksum(txdat = data.frame(x = log(dat$x + 1)), tydat = dat$y,
      bws = bw, bwtype = type, leave.one.out = TRUE)
    expect_equal(direct$ksum, reference$ksum, tolerance = 1e-14)
  }
})
