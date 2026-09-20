test_that("data-first native calls agree with named training roles", {
  withr::local_options(np.messages = FALSE)
  set.seed(4191)
  d <- data.frame(x = runif(24), y = rnorm(24), z = runif(24))
  d$y <- 2 * d$x + d$z + .3 * d$y
  cases <- list(
    npudens = list(tdat = d[c("x", "z")]),
    npudist = list(tdat = d[c("x", "z")]),
    npreg = list(txdat = d["x"], tydat = d$y),
    npcdens = list(txdat = d["x"], tydat = d["y"]),
    npcdist = list(txdat = d["x"], tydat = d["y"]),
    npindex = list(txdat = d[c("x", "z")], tydat = d$y),
    npplreg = list(txdat = d["x"], tydat = d$y, tzdat = d["z"]),
    npscoef = list(txdat = d["x"], tydat = d$y))
  bw.values <- function(b) if (inherits(b, "plbandwidth"))
    lapply(b$bw, function(component) component$bw) else b$bw
  for (family in names(cases)) {
    fit <- get(family, mode = "function")
    control <- if (family %in% c("npudens", "npudist"))
      list(bwmethod = "normal-reference") else
      list(nmulti = 1L, itmax = 1L, tol = .1, ftol = .1)
    if (family == "npscoef") control <- list(nmulti = 1L, optim.maxit = 10L, iterate = FALSE)
    control$se <- FALSE
    set.seed(4192)
    expected <- do.call(fit, c(cases[[family]], control))
    seed <- .Random.seed
    # All data-first defaults enter the retained bandwidth-object fit route.
    oracle <- fit(expected$bws, se = FALSE)
    for (wrapped in c(FALSE, TRUE)) {
      run <- if (wrapped) function(...) fit(...) else fit
      set.seed(4192)
      actual <- do.call(run, c(unname(cases[[family]]), control))
      expect_equal(fitted(actual), fitted(oracle), tolerance = 0,
                   info = paste(family, wrapped))
      expect_equal(bw.values(actual$bws), bw.values(expected$bws), tolerance = 0,
                   info = paste(family, wrapped))
      expect_identical(.Random.seed, seed, info = paste(family, wrapped))
      expect_equal(actual$bws[[".np.native.training"]],
                   expected$bws[[".np.native.training"]], tolerance = 0,
                   info = paste(family, wrapped))
    }
  }
})

test_that("unconditional data-first syntax accepts all native input shapes", {
  withr::local_options(np.messages = FALSE)
  x <- seq(-1, 1, length.out = 18)
  for (fit in list(npudens, npudist))
    for (value in list(x, matrix(x, ncol = 1L), data.frame(x = x))) {
      expected <- fit(tdat = value, bwmethod = "normal-reference")
      actual <- fit(value, bwmethod = "normal-reference")
      expect_identical(fitted(actual), fitted(expected))
      expect_identical(actual$bws$bw, expected$bws$bw)
      expect_identical(predict(actual, newdata = matrix(x[1:4], ncol = 1L)),
                       predict(expected, newdata = matrix(x[1:4], ncol = 1L)))
    }
})

test_that("automatic data-first promises are evaluated once", {
  withr::local_options(np.messages = FALSE)
  x <- data.frame(x = seq(-1, 1, length.out = 18))
  count <- 0L
  touch <- function() { count <<- count + 1L; x }
  wrapper <- function(...) npudens(...)
  fit <- wrapper(touch(), bwmethod = "normal-reference")
  expect_identical(count, 1L)
  expect_identical(fitted(fit), fitted(npudens(tdat = x, bwmethod = "normal-reference")))
})

test_that("optional smooth-coefficient z data do not displace manual bandwidths", {
  withr::local_options(np.messages = FALSE)
  x <- data.frame(x = seq(.1, 1, length.out = 18))
  y <- sin(x$x)
  a <- npscoef(.4, x, y, iterate = FALSE)
  b <- npscoef(bws = .4, txdat = x, tydat = y, iterate = FALSE)
  expect_identical(fitted(a), fitted(b))
  expect_identical(a$bws$bw, b$bws$bw)
})
