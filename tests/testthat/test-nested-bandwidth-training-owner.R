test_that("partially linear child bandwidths own their native samples", {
  withr::local_options(np.messages = FALSE)
  set.seed(1519)
  d <- data.frame(y = rnorm(24), x = rnorm(24), z = runif(24))
  check <- function(parent) {
    for (i in seq_along(parent$bw)) {
      child <- parent$bw[[i]]
      expect_null(environment(child$call))
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
  owner <- getFromNamespace(".np_bws_retain_formula_training", "np")
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
