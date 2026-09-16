test_that("copula formula entries share one constructor and fitting sample", {
  old <- options(np.messages = FALSE); on.exit(options(old), add = TRUE)
  set.seed(992L)
  d <- data.frame(x = runif(32), z = rnorm(32))
  count <- new.env(parent = emptyenv()); count$n <- 0L
  counted <- function(x) { count$n <- count$n + 1L; x }
  f <- ~ counted(x) + z
  dat <- d; names(dat)[1] <- "counted(x)"
  for (target in c("density", "distribution")) {
    constructor <- if (target == "density") npudensbw else npudistbw
    bw <- constructor(dat = dat, bws = c(.5, .6), bandwidth.compute = FALSE)
    for (grid in c(FALSE, TRUE)) {
      u <- if (grid) matrix(c(.2, .4, .6, .8), ncol = 2L) else NULL
      ref <- npcopula(bw, data = dat, u = u, se = TRUE, n.quasi.inv = 20)
      for (args in list(list(formula = f, bws = c(.5, .6), data = d),
          list(data = d, bws = c(.5, .6), formula = f))) {
        count$n <- 0L; rng <- .Random.seed
        value <- do.call(npcopula, c(args, list(target = target, u = u,
          evaluation = "sample", bandwidth.compute = FALSE, se = TRUE, n.quasi.inv = 20)))
        expect_identical(count$n, 1L)
        expect_identical(.Random.seed, rng)
        expect_identical(fitted(value), fitted(ref))
        expect_identical(se(value), se(ref))
        expect_identical(as.data.frame(value), as.data.frame(ref))
        expect_false(grepl("formula.state", paste(deparse(value$bws$call), collapse = "")))
        expect_false(grepl("getFromNamespace", paste(deparse(attr(value$bws$terms, "predvars")), collapse = "")))
      }
    }
  }
  count$n <- 0L
  positional <- npcopula(f, data = d, bws = c(.5, .6), target = "density",
    evaluation = "sample", bandwidth.compute = FALSE)
  expect_identical(count$n, 1L)
  bw <- npudensbw(dat = dat, bws = c(.5, .6), bandwidth.compute = FALSE)
  expect_identical(fitted(positional), fitted(npcopula(bw, data = dat)))
  count$n <- 0L
  value <- npcopula(f, data = d, target = "density", evaluation = "sample",
    bwmethod = "normal-reference")
  expect_identical(count$n, 1L)
  bw <- npudensbw(dat = dat, bwmethod = "normal-reference")
  expect_identical(fitted(value), fitted(npcopula(bw, data = dat)))
})

test_that("copula formula uses the first stochastic sample and subset", {
  old <- options(np.messages = FALSE); on.exit(options(old), add = TRUE)
  set.seed(993L)
  d <- data.frame(x = runif(32), z = rnorm(32))
  jittered <- function(x) x + rnorm(length(x), sd = .02)
  set.seed(994L)
  dat <- data.frame(jittered(d$x), d$z); names(dat) <- c("jittered(x)", "z")
  rng <- .Random.seed
  bw <- npudensbw(dat = dat, bws = c(.5, .6), bandwidth.compute = FALSE)
  ref <- npcopula(bw, data = dat, se = TRUE)
  set.seed(994L)
  value <- npcopula(bws = c(.5, .6), formula = ~ jittered(x) + z, data = d,
    target = "density", evaluation = "sample", bandwidth.compute = FALSE, se = TRUE)
  expect_identical(.Random.seed, rng)
  expect_identical(fitted(value), fitted(ref))
  expect_identical(se(value), se(ref))
  d$x[3] <- NA_real_
  for (target in c("density", "distribution")) {
    constructor <- if (target == "density") npudensbw else npudistbw
    bw <- constructor(~ x + z, data = d, subset = 2:25,
      bws = c(.5, .6), bandwidth.compute = FALSE)
    ref <- npcopula(bw, data = na.omit(d[2:25, ]))
    value <- npcopula(bws = c(.5, .6), formula = ~ x + z, data = d,
      subset = 2:25, target = target, evaluation = "sample", bandwidth.compute = FALSE)
    expect_identical(fitted(value), fitted(ref))
    expect_length(fitted(value), 23L)
  }
})

test_that("copula time-series formulas agree with independently indexed data", {
  old <- options(np.messages = FALSE); on.exit(options(old), add = TRUE)
  set.seed(995L)
  x <- ts(rnorm(36), frequency = 4)
  dat <- data.frame(as.numeric(x)[2:36], as.numeric(x)[1:35])
  names(dat) <- c("x", "lag(x, -1)")
  for (target in c("density", "distribution")) {
    constructor <- if (target == "density") npudensbw else npudistbw
    for (type in c("fixed", "generalized_nn", "adaptive_nn")) {
      widths <- if (type == "fixed") c(.8, .9) else c(15, 17)
      bw <- constructor(dat = dat, bws = widths, bwtype = type, bandwidth.compute = FALSE)
      ref <- npcopula(bw, data = dat)
      value <- npcopula(formula = ~ x + lag(x, -1), bws = widths, bwtype = type,
        target = target, evaluation = "sample", bandwidth.compute = FALSE)
      expect_identical(fitted(value), fitted(ref))
      expect_identical(as.data.frame(value), as.data.frame(ref))
      expect_length(fitted(value), 35L)
    }
  }
})

test_that("copula invalid formula input fails before bandwidth search", {
  old <- options(np.messages = FALSE); on.exit(options(old), add = TRUE)
  d <- data.frame(x = seq(.1, .9, length.out = 24), z = seq(.9, .1, length.out = 24), w = 1:24)
  for (target in c("density", "distribution")) {
    # A search-only invalid option is a negative control for early grid validation.
    expect_error(npcopula(~ x + z + w, data = d, target = target, bwmethod = "invalid"),
      "automatic copula probability grids")
    expect_error(npcopula(bws = c(.5, .6, .7), formula = ~ x + z + w,
      data = d, target = target, bwmethod = "invalid"), "automatic copula probability grids")
    expect_error(npcopula(~ absent + z, data = d, target = target), "absent")
    value <- npcopula(bws = c(.5, .6), formula = ~ x + z, data = d,
      target = target, evaluation = "sample", bandwidth.compute = FALSE)
    expect_length(fitted(value), 24L)
  }
})
