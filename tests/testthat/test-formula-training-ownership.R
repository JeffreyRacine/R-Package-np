test_that("formula bandwidths retain their own selected training values", {
  withr::local_options(np.messages = FALSE)
  set.seed(1031)
  d <- data.frame(y = rnorm(30), x = runif(30), z = rnorm(30),
                  keep = rep(c(TRUE, TRUE, FALSE), 10))
  d$x[5] <- NA_real_
  cases <- list(
    npreg = list(y ~ x, .4), npudens = list(~ x, .4),
    npudist = list(~ x, .4), npcdens = list(y ~ x, c(.4, .5)),
    npcdist = list(y ~ x, c(.4, .5)),
    npplreg = list(y ~ x | z, matrix(.5, 2, 1)),
    npscoef = list(y ~ x | z, .5), npindex = list(y ~ x + z, c(1, .5, .4)))
  for (family in names(cases)) {
    constructor <- get(paste0(family, "bw"))
    fit <- get(family)
    spec <- cases[[family]]
    dat <- d
    policy <- na.exclude
    bw <- constructor(spec[[1]], data = dat, subset = keep, na.action = policy,
                      bws = spec[[2]], bandwidth.compute = FALSE)
    ref <- fitted(fit(bw, se = FALSE))
    dat <- transform(d, y = y + 17, x = x + 4, z = z - 8)
    policy <- function(x) stop("rebound policy must not be called")
    expect_identical(fitted(fit(bw, se = FALSE)), ref, info = family)
    restored <- unserialize(serialize(bw, NULL))
    expect_identical(fitted(fit(restored, se = FALSE)), ref, info = family)
    expect_identical(bw[[".np.formula.training"]]$na.action, stats::na.exclude)
    # Explicit training data replacement remains supported, using the
    # retained NA policy rather than the rebound constructor variable.
    control <- constructor(spec[[1]], data = dat, subset = keep,
      na.action = na.exclude, bws = spec[[2]], bandwidth.compute = FALSE)
    expect_identical(fitted(fit(bw, data = dat, se = FALSE)),
                     fitted(fit(control, se = FALSE)), info = family)
  }
})

test_that("retained frames do not re-run transforms, data or NA factories", {
  withr::local_options(np.messages = FALSE)
  set.seed(1032)
  d <- data.frame(x = runif(24), y = rnorm(24))
  d$y[4] <- NA_real_
  counts <- c(data = 0L, policy = 0L, action = 0L, transform = 0L)
  values <- function() { counts["data"] <<- counts["data"] + 1L; d }
  policy <- function() {
    counts["policy"] <<- counts["policy"] + 1L
    function(x) { counts["action"] <<- counts["action"] + 1L; na.exclude(x) }
  }
  jitter.once <- function(x) {
    counts["transform"] <<- counts["transform"] + 1L
    x + runif(length(x)) / 10
  }
  bw <- npregbw(y ~ jitter.once(x), data = values(), na.action = policy(),
                bws = .4, bandwidth.compute = FALSE)
  expect_identical(unname(counts), rep(1L, 4))
  seed <- .Random.seed
  ref <- npreg(bw)
  again <- npreg(bw)
  expect_identical(fitted(again), fitted(ref))
  expect_identical(unname(counts), rep(1L, 4))
  expect_identical(.Random.seed, seed)
  x <- d$x
  y <- d$y
  b0 <- npregbw(y ~ x, bws = .4, bandwidth.compute = FALSE)
  a <- fitted(npreg(b0))
  x <- x + 8
  y <- y + 13
  expect_identical(fitted(npreg(b0)), a)
  # Explicit data environments are materialized, not retained by reference.
  e <- list2env(d, parent = baseenv())
  be <- npregbw(y ~ x, data = e, bws = .4, bandwidth.compute = FALSE)
  a <- fitted(npreg(be))
  e$y <- rep(11, 24)
  expect_identical(fitted(npreg(be)), a)
})

test_that("single-index plotting uses the same saved training owner", {
  withr::local_options(np.messages = FALSE)
  set.seed(1033)
  d <- data.frame(x = runif(30), z = rnorm(30), y = rnorm(30))
  d$y[4] <- NA_real_
  wrapper <- function(...) npindexbw(...)
  b <- wrapper(y ~ x + z, data = d, na.action = na.exclude,
               bws = c(1, .2, .5), bandwidth.compute = FALSE)
  ref <- plot(b, output = "data", errors = "none", neval = 3L)
  d$y <- rep(11, 30)
  a <- plot(b, output = "data", errors = "none", neval = 3L)
  expect_identical(a, ref)
})

test_that("character NA policies use the model-frame owner, not local names", {
  withr::local_options(np.messages = FALSE, na.action = "na.omit")
  d <- data.frame(y = sin(1:20), x = seq_len(20) / 20)
  d$y[c(3, 7)] <- NA_real_
  hits <- 0L
  na.omit <- na.exclude <- function(x, ...) { hits <<- hits + 1L; x[1:4, ] }
  for (policy in c("na.omit", "na.exclude")) {
    b <- npregbw(y ~ x, data = d, na.action = policy, bws = .4, bandwidth.compute = FALSE)
    oracle <- stats::model.frame(y ~ x, data = d, na.action = policy)
    expect_identical(b[[".np.formula.training"]]$frame, oracle)
  }
  b <- npregbw(y ~ x, data = d, bws = .4, bandwidth.compute = FALSE)
  expect_equal(b$nobs, 18)
  expect_identical(hits, 0L)
  # A function supplied intentionally still belongs to the user.
  custom <- function(x) { hits <<- hits + 1L; stats::na.exclude(x) }
  b <- npregbw(y ~ x, data = d, na.action = custom, bws = .4, bandwidth.compute = FALSE)
  expect_identical(b[[".np.formula.training"]]$na.action, custom)
  expect_identical(hits, 1L)
  npreg(b)
  expect_identical(hits, 1L)
})

test_that("replacement-data policy overrides are honored and retained across families", {
  withr::local_options(np.messages = FALSE)
  set.seed(19312)
  d <- data.frame(y = rnorm(30), x = runif(30), z = rnorm(30))
  newer <- transform(d, y = y + 2)
  newer$x[c(4, 13)] <- NA_real_
  cases <- list(npreg = list(y ~ x, .5), npudens = list(~ x, .5),
    npudist = list(~ x, .5), npcdens = list(y ~ x, c(.5, .7)),
    npcdist = list(y ~ x, c(.5, .7)), npindex = list(y ~ x + z, c(1, .3, .5)),
    npplreg = list(y ~ x | z, matrix(.6, 2, 1)), npscoef = list(y ~ x | z, .6))
  for (family in names(cases)) {
    ctor <- get(paste0(family, "bw")); fit <- get(family); spec <- cases[[family]]
    b <- ctor(spec[[1L]], data = d, bws = spec[[2L]], bandwidth.compute = FALSE)
    expect_error(fit(b, data = newer, na.action = stats::na.fail, se = FALSE),
                 "missing values", info = family)
    for (policy in list(stats::na.exclude, stats::na.omit, stats::na.pass, NULL)) {
      expected.b <- ctor(spec[[1L]], data = newer, na.action = policy,
                         bws = spec[[2L]], bandwidth.compute = FALSE)
      a <- fit(b, data = newer, na.action = policy, se = FALSE)
      expected <- fit(expected.b, se = FALSE)
      expect_equal(fitted(a), fitted(expected), tolerance = 1e-12, info = family)
      expect_identical(a$bws[[".np.formula.training"]]$na.action, policy, info = family)
      expect_equal(fitted(fit(unserialize(serialize(a$bws, NULL)), se = FALSE)),
                   fitted(expected), tolerance = 1e-12, info = family)
      expect_identical(nrow(a$bws[[".np.formula.training"]]$frame), 28L, info = family)
      expect_identical(as.integer(attr(a$bws[[".np.formula.training"]]$frame, "na.action")),
                       c(4L, 13L), info = family)
    }
  }
})

test_that("training completion records native exclusions without replaying a policy", {
  withr::local_options(np.messages = FALSE)
  d <- data.frame(x = seq_len(20) / 20, y = sin(1:20))
  d$x[c(2, 8)] <- NA_real_
  calls <- 0L
  policy <- function(x) { calls <<- calls + 1L; x }
  a <- npreg(y ~ x, data = d, na.action = policy, bws = .4, bandwidth.compute = FALSE)
  expect_identical(calls, 1L)
  expect_equal(a$rows.omit, c(2, 8))
  expect_equal(fitted(a), fitted(npreg(y ~ x, data = stats::na.omit(d),
    bws = .4, bandwidth.compute = FALSE)))
  expect_equal(fitted(npreg(a$bws)), fitted(a))
  expect_identical(calls, 1L)
})

test_that("derived conditional and regression consumers honor replacement NA policy", {
  withr::local_options(np.messages = FALSE)
  d <- data.frame(x = seq_len(24) / 24, y = sin(1:24), group = factor(rep(1:2, 12)))
  r <- npregbw(y ~ x, data = d, bws = .4, bandwidth.compute = FALSE)
  q <- npcdistbw(y ~ x, data = d, bws = c(.4, .6), bandwidth.compute = FALSE)
  c <- npcdensbw(group ~ x, data = d, bws = c(.4, .2), bandwidth.compute = FALSE)
  d$x[4] <- NA_real_
  expect_error(npreghat(r, data = d, na.action = stats::na.fail), "missing values")
  expect_error(npsigtest(r, data = d, na.action = stats::na.fail, B = 9L), "missing values")
  expect_error(npqreg(q, data = d, na.action = stats::na.fail), "missing values")
  expect_error(npconmode(c, data = d, na.action = stats::na.fail), "missing values")
})
