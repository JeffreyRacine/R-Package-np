test_that("raw index formulas consume NA policy before native handoff", {
  old <- options(np.messages = FALSE, na.action = "na.omit")
  on.exit(options(old), add = TRUE)
  set.seed(423)
  d <- data.frame(x = runif(70), z = runif(70))
  d$y <- sin(3 * (d$x + .5 * d$z)) + rnorm(70, sd = .2)
  d$x[c(3, 9)] <- NA
  nw <- function(train, ev = train) {
    train <- train[complete.cases(train), ]
    v <- train$x + .5 * train$z
    vapply(ev$x + .5 * ev$z, function(u) {
      k <- dnorm((u - v) / .2)
      sum(k * train$y) / sum(k)
    }, 0)
  }
  b <- c(1, .5, .2)
  for (policy in list(na.omit, "na.omit", na.exclude, NULL, na.pass)) {
    f <- npindex(y ~ x + z, data = d, bws = b, se = FALSE, na.action = policy)
    exclude <- identical(policy, na.exclude)
    expected <- nw(d, if (exclude) d else d[complete.cases(d), ])
    expect_equal(fitted(f), expected, tolerance = 1e-12, ignore_attr = TRUE)
    expect_equal(f$ntrain, 68L)
  }
  expect_error(npindex(y ~ x + z, data = d, bws = b, se = FALSE,
                       na.action = na.fail), "missing values")
  f <- npindex(y ~ x + z, data = d, subset = z > .25, bws = b,
               na.action = na.omit, se = FALSE)
  tr <- subset(d, z > .25); tr <- tr[complete.cases(tr), ]
  expect_equal(fitted(f), nw(tr), tolerance = 1e-12, ignore_attr = TRUE)
  e <- data.frame(x = c(.1, NA, .5, .8), z = c(.2, .3, .7, .9),
                  y = c(0, 1, NA, 2))
  f <- npindex(y ~ x + z, data = d, newdata = e, bws = b,
               na.action = na.exclude, se = FALSE)
  expect_equal(fitted(f), nw(d, e), tolerance = 1e-12, ignore_attr = TRUE)
  expect_equal(f$ntrain, 68L)
  expect_equal(f$rows.omit, 2L)
  expect_true(is.finite(fitted(f)[3]))
  expect_error(npindex(y ~ x + z, data = d[complete.cases(d), ],
                       newdata = e, bws = b, se = FALSE, na.action = na.fail),
               "missing values")
  calls <- new.env(); calls$force <- calls$apply <- calls$data <- 0L
  policy <- function(x) { calls$apply <- calls$apply + 1L; na.exclude(x) }
  f <- npindex(y ~ x + z,
    data = { calls$data <- calls$data + 1L; d }, newdata = e, bws = b,
    na.action = { calls$force <- calls$force + 1L; policy }, se = FALSE)
  expect_equal(c(calls$data, calls$force, calls$apply), c(1L, 1L, 2L))
  expect_equal(fitted(f), nw(d, e), tolerance = 1e-12, ignore_attr = TRUE)
  for (args in list(list(y ~ x + z), list(formula = y ~ x + z))) {
    f <- do.call(npindex, c(args, list(data = d, bws = b,
                                     na.action = na.omit, se = FALSE)))
    expect_equal(fitted(f), nw(d, d[complete.cases(d), ]),
                 tolerance = 1e-12, ignore_attr = TRUE)
  }
  # Direct native calls still refuse a formula-only policy.
  expect_error(npindexbw(xdat = d[c("x", "z")], ydat = d$y,
                         bws = b, bandwidth.compute = FALSE, na.action = na.fail),
               "supported only with the formula interface", fixed = TRUE)
})

test_that("raw index formulas preserve binary models and derivative routes", {
  old <- options(np.messages = FALSE); on.exit(options(old), add = TRUE)
  set.seed(152); d <- data.frame(x = rnorm(85), z = rnorm(85))
  d$y <- sin(d$x + .35*d$z) + rnorm(85, sd = .3)
  d$x[4] <- NA
  for (method in c("ichimura", "kleinspady")) {
    dd <- d
    if (method == "kleinspady") dd$y <- as.integer(dd$y > 0)
    b <- npindexbw(y ~ x + z, data = dd, bws = c(1, .35, .7),
                   bandwidth.compute = FALSE, method = method, na.action = na.exclude)
    ref <- npindex(b, gradients = TRUE, se = TRUE)
    fit <- npindex(y ~ x + z, data = dd, bws = c(1, .35, .7),
                   method = method, na.action = na.exclude, gradients = TRUE, se = TRUE)
    for (key in c("mean", "merr", "grad", "gerr", "ntrain", "rows.omit"))
      expect_equal(fit[[key]], ref[[key]], info = key)
  }
})
