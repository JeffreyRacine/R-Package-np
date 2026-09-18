test_that("partially linear child bandwidths own their native samples", {
  withr::local_options(np.messages = FALSE)
  set.seed(1519)
  d <- data.frame(y = rnorm(24), x = rnorm(24), z = runif(24))
  check <- function(parent) {
    for (i in seq_along(parent$bw)) {
      child <- parent$bw[[i]]
      if (!is.null(child$call)) expect_null(environment(child$call))
      values <- child[[".np.native.training", exact = TRUE]]
      expect_named(values, c("xdat", "ydat"))
      expect_equal(values$xdat, d["z"], ignore_attr = TRUE)
      expect_equal(unname(values$ydat), if (i == 1L) d$y else d$x)
      control <- npreg(child, txdat = d["z"],
                        tydat = if (i == 1L) d$y else d$x, se = TRUE)
      restored <- unserialize(serialize(child, NULL))
      actual <- npreg(restored, se = TRUE)
      expect_identical(fitted(actual), fitted(control))
      expect_identical(se(actual), se(control))
      expect_identical(predict(actual, exdat = d[1:3, "z", drop = FALSE]),
                       predict(control, exdat = d[1:3, "z", drop = FALSE]))
      # A caller with the old internal names must not replace retained data.
      zdat <- data.frame(z = rep(80, 24))
      ydat <- rep(-9, 24)
      expect_identical(fitted(npreg(restored)), fitted(control))
    }
  }
  b <- npplregbw(y ~ x | z, data = d, bws = matrix(.5, 2, 1),
                  bandwidth.compute = FALSE)
  check(b)
  # Exercise the ordinary reselection owner without changing its arithmetic.
  searched <- npplregbw(xdat = d["x"], ydat = d$y, zdat = d["z"],
                         bws = b, nmulti = 1L, itmax = 2L)
  check(searched)
  # Different fixed child degrees exercise the shared builder beyond LC.
  polynomial <- npplregbw(y ~ x | z, data = d, bws = matrix(.5, 2, 1),
    regtype = "lp", degree = list(1L, 2L), bandwidth.compute = FALSE)
  check(polynomial)
  a <- plot(b$bw[[1L]], errors = "none", output = "data", neval = 3L)
  zdat <- data.frame(z = rep(80, 24)); ydat <- rep(-9, 24)
  expect_identical(plot(b$bw[[1L]], errors = "none",
                       output = "data", neval = 3L), a)
})

test_that("local-smoothing quantile children retain distinct pilot responses", {
  withr::local_options(np.messages = FALSE)
  set.seed(17653)
  x <- data.frame(x = seq(-1, 1, length.out = 24))
  y <- sin(x$x) + rnorm(24)
  model <- nplsqreg(txdat = x, tydat = y, bws = .7, delta = .5,
    tau = c(.25, .75), bandwidth.compute = FALSE, nomad = FALSE,
    regtype = "ll", se = TRUE)
  for (one in model$tau.fits) {
    expect_true(is.symbol(one$call[[1L]]))
    expect_null(environment(one$call))
    expect_true(is.symbol(one$bws$call[[1L]]))
    expect_null(environment(one$bws$call))
    expect_identical(one$fit$bws[[".np.native.training"]]$ydat, one$bws$qdat)
    children <- list(one$fit)
    if (!is.null(one$bws$mean.fit)) {
      expect_identical(one$bws$mean.fit$bws[[".np.native.training"]]$ydat, y)
      expect_identical(one$bws$scale.fit$bws[[".np.native.training"]]$ydat,
                       (y - as.numeric(fitted(one$bws$mean.fit)))^2)
      children <- c(children, list(one$bws$mean.fit, one$bws$scale.fit))
    }
    for (child in children) {
      expect_null(environment(child$call))
      if (!is.null(child$bws$call)) expect_null(environment(child$bws$call))
      restored <- unserialize(serialize(child$bws, NULL))
      expect_equal(fitted(npreg(restored)), fitted(child), tolerance = 1e-14)
      # Internal native bandwidths without an original call cannot be updated
      # by evaluating package-invented symbols in the user's frame.
      if (is.null(child$bws$call)) {
        xdat <- data.frame(x = seq(1000, 2000, length.out = 24))
        ydat <- seq_len(24)
        expect_error(update(restored), "need an object with call component")
        expect_equal(fitted(npreg(restored)), fitted(child), tolerance = 1e-14)
        values <- restored[[".np.native.training"]]
        expect_equal(as.numeric(npreghat(restored)),
                     as.numeric(npreghat(restored, txdat = values$xdat)),
                     tolerance = 1e-14)
        expect_equal(plot(restored, errors = "none", output = "data", neval = 3L),
          plot(restored, xdat = values$xdat, ydat = values$ydat,
               errors = "none", output = "data", neval = 3L))
        set.seed(7301)
        actual <- npsigtest(restored, B = 9L)
        set.seed(7301)
        expected <- npsigtest(restored, xdat = values$xdat, ydat = values$ydat, B = 9L)
        expect_identical(actual$P, expected$P)
      }
    }
  }
})

test_that("retaining a call-less native child never fabricates executable provenance", {
  retain <- getFromNamespace(".np_bws_retain_native_training", "npRmpi")
  xdat <- data.frame(x = seq_len(6))
  ydat <- 6:1
  b <- retain(list(), xdat = xdat, ydat = ydat)
  expect_null(b$call)
  expect_identical(b[[".np.native.training"]], list(xdat = xdat, ydat = ydat))
  expect_error(update(b), "need an object with call component")
})

test_that("child sample retention does not strip user formula environments", {
  withr::local_options(np.messages = FALSE)
  owner <- new.env(parent = globalenv())
  owner$shift <- function(x) x + .1
  owner$d <- data.frame(y = seq_len(16), x = seq_len(16) / 16)
  f <- y ~ shift(x)
  environment(f) <- owner
  b <- npregbw(f, data = owner$d, bws = .5, bandwidth.compute = FALSE)
  expect_identical(environment(b$formula), owner)
  expect_null(b[[".np.native.training", exact = TRUE]])
  a <- npreg(b)
  expect_equal(predict(a, newdata = owner$d[1:3, , drop = FALSE]),
               fitted(npreg(b, newdata = owner$d[1:3, , drop = FALSE])))
})

test_that("formula sample retention releases only internal call owners", {
  owner <- getFromNamespace(".np_bws_retain_formula_training", "npRmpi")
  user <- new.env(parent = globalenv())
  formula <- y ~ x; environment(formula) <- user
  frame <- data.frame(y = 1:4, x = 4:1)
  b <- list(call = quote(npregbw(y ~ x)), formula = formula)
  internal <- new.env(parent = environment(owner))
  for (call.owner in list(user, internal, emptyenv())) {
    environment(b$call) <- call.owner
    out <- owner(b, frame, stats::na.exclude)
    expect_identical(environment(out$call),
      if (identical(call.owner, internal)) user else call.owner)
    expect_identical(environment(out$formula), user)
    expect_identical(out[[".np.formula.training"]]$frame, frame)
  }
})
