test_that("canonical frames evaluate expressions once and preserve prediction metadata", {
  frame <- getFromNamespace(".np_formula_model_frame", "npRmpi")
  set.seed(842L)
  d <- data.frame(x = rnorm(40), y = rnorm(40), z = rnorm(40),
                  f = factor(rep(letters[1:4], 10)))
  count <- 0L
  counted <- function(x) {count <<- count + 1L; x}
  actual <- frame(y ~ counted(x) + z + f, data = d)
  expect_identical(count, 1L)
  expect_identical(actual, model.frame(y ~ counted(x) + z + f, data = d))
  for (formula in list(y ~ poly(x, 2) + f, y ~ splines::ns(x, df = 3) + f,
                       y ~ ., y ~ log(abs(x) + 1) + f, y ~ I(cbind(x, z)))) {
    actual <- frame(formula, data = d)
    oracle <- model.frame(formula, data = d)
    expect_identical(actual, oracle)
    new <- d[4:9, , drop = FALSE]
    expect_identical(frame(attr(actual, "terms"), data = new),
                     model.frame(attr(oracle, "terms"), data = new))
  }
  keep <- FALSE
  d$keep <- rep(c(TRUE, FALSE), 20)
  expect_identical(frame(y ~ x + f, data = d, subset = keep),
                   model.frame(y ~ x + f, data = d, subset = keep))
  d$x[3] <- NA
  expect_identical(frame(y ~ x + f, data = d, na.action = na.exclude),
                   model.frame(y ~ x + f, data = d, na.action = na.exclude))
  tt <- terms(y ~ x)
  attr(tt, "predvars") <- quote(list(y, x + 1))
  expect_identical(frame(tt, data = d), model.frame(tt, data = d))
  # Only the exact historical package-owned wrapper may be removed.
  for (origin in c("np", "npRmpi")) {
    tt <- terms(y ~ counted(x))
    attr(tt, "predvars") <- substitute(
      utils::getFromNamespace(".np_formula_align_values", ORIGIN)(list(y, counted(x))),
      list(ORIGIN = origin))
    count <- 0L
    actual <- frame(tt, data = d)
    expect_identical(count, 1L)
    expect_identical(attr(attr(actual, "terms"), "predvars"), quote(list(y, counted(x))))
  }
  y <- ts(1:12, start = c(2000, 1), frequency = 4)
  x <- ts(cbind(a = 101:112, b = 201:212), start = c(2000, 1), frequency = 4)
  actual <- frame(y ~ lag(x, -1) + lag(y, -1))
  expect_identical(dim(actual[[2L]]), c(11L, 2L))
  expect_identical(as.numeric(actual[[3L]]), as.numeric(y)[1:11])
  expect_identical(attr(attr(actual, "terms"), "predvars"),
                   quote(list(y, lag(x, -1), lag(y, -1))))
})

test_that("one regression transaction reuses its sample and independent refits do not", {
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  set.seed(843L)
  d <- data.frame(x = rnorm(40), y = rnorm(40))
  calls <- c(x = 0L, y = 0L)
  counted.x <- function(x) {calls["x"] <<- calls["x"] + 1L; x}
  counted.y <- function(y) {calls["y"] <<- calls["y"] + 1L; y}
  f <- counted.y(y) ~ counted.x(x)
  one <- npreg(f, data = d, bws = .5)
  expect_identical(unname(calls), c(1L, 1L))
  calls[] <- 0L
  two <- npreg(formula = f, data = d, bws = .5)
  expect_identical(unname(calls), c(1L, 1L))
  expect_identical(one$mean, two$mean)
  calls[] <- 0L
  separate <- npreg(bws = one$bws)
  expect_identical(unname(calls), c(1L, 1L))
  expect_identical(one$mean, separate$mean)
  calls[] <- 0L
  auto <- npreg(f, data = d, nmulti = 1L)
  expect_identical(unname(calls), c(1L, 1L))
  expect_false(".np.formula.state" %in% names(auto$bws$call))
  expect_false(".np.formula.state" %in% names(auto$bws$call$...))
  expect_false(".np.formula.state" %in% names(auto$call$...))
  expect_identical(attr(one$bws$terms, "predvars"),
                   quote(list(counted.y(y), counted.x(x))))
  # A stochastic formula must use the same draw for constructor and fit.
  jittered <- function(x) x + runif(length(x), -.1, .1)
  set.seed(844L)
  actual <- npreg(y ~ jittered(x), data = d, bws = .5)
  rng <- .Random.seed
  set.seed(844L)
  x.oracle <- data.frame(x = jittered(d$x))
  expected <- npreg(txdat = x.oracle, tydat = d$y, bws = .5)
  expect_identical(.Random.seed, rng)
  expect_equal(actual$mean, expected$mean, tolerance = 1e-14)
})

test_that("formula preparation contexts are single-use and never global caches", {
  store <- getFromNamespace(".np_formula_frame_store", "npRmpi")
  take <- getFromNamespace(".np_formula_frame_take", "npRmpi")
  state <- new.env(parent = emptyenv())
  value <- data.frame(x = 1:3)
  store(state, value)
  expect_error(store(state, value), "invalid formula preparation context")
  expect_identical(take(state), value)
  expect_identical(ls(state, all.names = TRUE), character())
  expect_error(take(state), "no training frame")
  expect_error(store(new.env(), value), "invalid formula preparation context")
})

test_that("legacy alignment wrappers are normalized before response deletion", {
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  set.seed(845L)
  d <- data.frame(x = rnorm(40), y = rnorm(40))
  bw <- npregbw(y ~ x, data = d, bws = .5, bandwidth.compute = FALSE)
  new <- d[1:4, , drop = FALSE]
  expected <- npreg(bws = bw, newdata = new)$mean
  for (origin in c("np", "npRmpi")) {
    legacy <- bw
    attr(legacy$terms, "predvars") <- substitute(
      utils::getFromNamespace(".np_formula_align_values", ORIGIN)(list(y, x)),
      list(ORIGIN = origin))
    expect_equal(npreg(bws = legacy, newdata = new)$mean, expected, tolerance = 1e-14)
    expect_equal(as.vector(npreghat(bws = legacy, newdata = new, output = "apply")),
                 expected, tolerance = 1e-14)
  }
})
