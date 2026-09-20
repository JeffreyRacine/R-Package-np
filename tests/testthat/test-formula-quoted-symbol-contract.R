test_that("formula labels resolve symbol names without evaluating expressions", {
  decode <- getFromNamespace(".np_formula_term_names", "np")
  expect_identical(decode(c("x", "`annual income`", "I(x + 1)", "`a`:`b`")),
    c("x", "annual income", "I(x + 1)", "`a`:`b`"))
  split <- getFromNamespace("explodePipe", "np")
  expect_identical(split(y ~ I(x + 1) + `x+y` | `a|b`),
    list("y", c("I(x + 1)", "`x+y`"), "`a|b`"))
})

test_that("quoted regression names survive fitting, refitting, prediction and tests", {
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  set.seed(427); d <- data.frame(x = runif(24), y = rnorm(24))
  control <- npreg(y ~ x, data = d, bws = .4)
  for (name in c("annual income", "if", "a`b", "x+y")) {
    dat <- d; names(dat)[1L] <- name
    f <- as.call(list(quote(`~`), quote(y), as.name(name)))
    class(f) <- "formula"; environment(f) <- environment()
    fit <- npreg(formula = f, data = dat, bws = .4)
    expect_identical(fitted(fit), fitted(control))
    expect_identical(fitted(npreg(fit$bws)), fitted(control))
    expect_identical(predict(fit, newdata = dat[1:4, ]),
                     predict(control, newdata = d[1:4, ]))
    expect_identical(fit$bws$xnames, name)
    expect_equal(unname(npsigtest(fit, B = 9)$P),
                 unname(npsigtest(control, B = 9)$P), tolerance = 2e-12)
  }
})

test_that("pipe-role labels preserve quoted punctuation and transformations", {
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  set.seed(428); d <- data.frame(x = runif(24), y = rnorm(24), z = rnorm(24))
  named <- d; names(named) <- c("x+y", "out come", "a|b")
  for (family in c("npplreg", "npscoef")) {
    fun <- get(family, asNamespace("np"))
    bw <- if (family == "npplreg") matrix(.4, 2, 1) else .4
    a <- do.call(fun, list(formula = y ~ x | z, data = d, bws = bw, bandwidth.compute = FALSE))
    b <- do.call(fun, list(formula = `out come` ~ `x+y` | `a|b`,
      data = named, bws = bw, bandwidth.compute = FALSE))
    expect_identical(fitted(a), fitted(b))
    expect_identical(unname(predict(a, newdata = d[1:4, ])),
                     unname(predict(b, newdata = named[1:4, ])))
  }
})
