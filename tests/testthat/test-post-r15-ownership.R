test_that("replacement fits and native bandwidths retain the actual training sample", {
  withr::local_options(np.messages = FALSE)
  set.seed(19310)
  d <- data.frame(x = runif(36, -1, 1), z = rnorm(36), y = rnorm(36))
  newer <- transform(d, x = x + .2, z = z / 2, y = 3 + x + .2 * y)
  cases <- list(npreg = list(y ~ x, .5), npudens = list(~ x, .5),
    npudist = list(~ x, .5), npcdens = list(y ~ x, c(.5, .7)),
    npcdist = list(y ~ x, c(.5, .7)), npindex = list(y ~ x + z, c(1, .3, .5)),
    npplreg = list(y ~ x | z, matrix(.6, 2, 1)), npscoef = list(y ~ x | z, .6))
  native <- function(family, data) {
    if (family %in% c("npudens", "npudist")) return(list(dat = data["x"]))
    out <- list(xdat = data[if (family == "npindex") c("x", "z") else "x"], ydat = data$y)
    if (family %in% c("npplreg", "npscoef")) out$zdat <- data["z"]
    out
  }
  for (family in names(cases)) {
    bwfun <- get(paste0(family, "bw"))
    fitfun <- get(family)
    spec <- cases[[family]]
    bw <- bwfun(spec[[1]], data = d, bws = spec[[2]], bandwidth.compute = FALSE)
    before <- bw
    fresh <- bwfun(spec[[1]], data = newer, bws = spec[[2]], bandwidth.compute = FALSE)
    fit <- fitfun(bw, data = newer, se = FALSE)
    oracle <- fitfun(fresh, se = FALSE)
    expect_equal(fitted(fit), fitted(oracle), tolerance = 1e-12, info = family)
    expect_equal(predict(fit, newdata = newer[1:5, ]), predict(oracle, newdata = newer[1:5, ]),
      tolerance = 1e-12, info = family)
    expect_equal(fitted(fitfun(fit$bws, se = FALSE)), fitted(oracle),
      tolerance = 1e-12, info = family)
    expect_identical(bw, before, info = family)
    expect_null(fit$bws[[".np.native.training"]])
    training <- native(family, d)
    b <- do.call(bwfun, c(training, list(bws = spec[[2]], bandwidth.compute = FALSE)))
    control <- fitted(fitfun(b, se = FALSE))
    training <- native(family, newer)
    restored <- b
    environment(restored$call) <- NULL
    expect_identical(fitted(fitfun(unserialize(serialize(restored, NULL)), se = FALSE)), control)
    names(training) <- c(dat = "tdat", xdat = "txdat", ydat = "tydat", zdat = "tzdat")[names(training)]
    explicit <- do.call(fitfun, c(list(bws = bw, se = FALSE), training))
    expect_equal(predict(explicit, newdata = newer[1:5, ]), predict(oracle, newdata = newer[1:5, ]),
      tolerance = 1e-12, info = family)
    named <- do.call(fitfun, c(list(bws = spec[[2]], bandwidth.compute = FALSE, se = FALSE), training))
    positional <- do.call(fitfun, c(list(spec[[2]]), training,
                                    list(bandwidth.compute = FALSE, se = FALSE)))
    expect_equal(fitted(positional), fitted(named), tolerance = 0, info = family)
    bw.values <- function(b) if (inherits(b, "plbandwidth"))
      lapply(b$bw, function(component) component$bw) else b$bw
    expect_equal(bw.values(positional$bws), bw.values(named$bws), tolerance = 0, info = family)
  }
})

test_that("derived refits use their own responses and not a competing formula snapshot", {
  withr::local_options(np.messages = FALSE)
  set.seed(19311)
  d <- data.frame(x = runif(32, -1, 1), y = rnorm(32))
  newer <- transform(d, y = 4 + x + .2 * y)
  b <- npcdistbw(y ~ x, data = d, bws = c(.5, .7), bandwidth.compute = FALSE)
  ref <- npcdistbw(y ~ x, data = newer, bws = c(.5, .7), bandwidth.compute = FALSE)
  a <- npqreg(b, data = newer, tau = c(.25, .75))
  r <- npqreg(ref, tau = c(.25, .75))
  expect_equal(predict(a, newdata = newer[1:4, ]), predict(r, newdata = newer[1:4, ]), tolerance = 1e-12)
  pb <- npregbw(y ~ x, data = d, bws = .5, regtype = "ll", bandwidth.compute = FALSE)
  q <- nplsqreg(txdat = d["x"], tydat = d$y, bws = .5, delta = .5,
    bandwidth.compute = FALSE, regtype = "ll", pilot.args = list(bws = pb))
  for (name in c("mean.fit", "scale.fit")) {
    child <- q$bws[[name]]
    expect_equal(fitted(npreg(child$bws)), fitted(child), tolerance = 1e-12)
    expect_true(xor(is.null(child$bws[[".np.native.training"]]),
                    is.null(child$bws[[".np.formula.training"]])))
  }
})

test_that("native refits do not reevaluate captured promises and formula shorthand works", {
  withr::local_options(np.messages = FALSE)
  d <- data.frame(x = seq(-1, 1, length.out = 24), y = sin(seq_len(24)))
  nx <- ny <- 0L
  b <- npregbw(xdat = { nx <- nx + 1L; d["x"] },
                ydat = { ny <- ny + 1L; d$y }, bws = .4, bandwidth.compute = FALSE)
  # Construction dispatch has a separately recorded incumbent promise issue.
  # This retained-data contract prohibits any additional evaluations on refit.
  forced <- c(nx, ny)
  a <- fitted(npreg(b))
  d$y <- d$y + 10
  expect_identical(fitted(npreg(b)), a)
  expect_identical(c(nx, ny), forced)
  wrapper <- function(...) npreg(...)
  direct <- npreg(y ~ x, data = d, bws = .4, bandwidth.compute = FALSE)
  expect_identical(fitted(wrapper(y ~ x, data = d, bws = .4, bandwidth.compute = FALSE)),
                   fitted(direct))
})
