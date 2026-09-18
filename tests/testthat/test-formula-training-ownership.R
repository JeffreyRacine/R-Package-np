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
