test_that("derived formula transactions preserve stochastic inputs and inference", {
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  set.seed(891L)
  d <- data.frame(y = rnorm(24), cy = factor(rep(1:2, 12)), x = runif(24))
  jittered <- function(x) x + runif(length(x), -.02, .02)
  for (family in c("npqreg", "npconmode")) {
    fun <- get(family)
    f <- if (family == "npqreg") y ~ jittered(x) else cy ~ jittered(x)
    h <- if (family == "npqreg") c(.7, .6) else c(.2, .6)
    controls <- if (family == "npqreg") list(tau = c(.25, .75)) else
      list(probabilities = TRUE, level = "1")
    set.seed(892L)
    actual <- do.call(fun, c(list(formula = f, data = d, bws = h,
      se = TRUE, gradients = TRUE), controls))
    rng <- .Random.seed
    set.seed(892L)
    x <- data.frame(x = jittered(d$x))
    response <- if (family == "npqreg") d$y else d$cy
    expected <- do.call(fun, c(list(bws = actual$bws, txdat = x, tydat = response,
      se = TRUE, gradients = TRUE), controls))
    expect_identical(.Random.seed, rng)
    fields <- if (family == "npqreg") c("quantile", "quanterr", "quantgrad", "quantgerr") else
      c("conmode", "probabilities", "probability.gradients", "probability.errors")
    for (field in fields) expect_equal(actual[[field]], expected[[field]], tolerance = 1e-12)
    expect_true(all(fields %in% names(actual)))
    expect_true(all(vapply(actual[fields], length, integer(1L)) > 0L))
  }
})

test_that("derived one-call frames retain subset and evaluation roles", {
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  set.seed(893L)
  d <- data.frame(y = rnorm(30), cy = factor(rep(1:2, 15)), x = runif(30))
  d$x[4] <- NA_real_
  keep <- seq_len(30) != 9L
  clean <- na.omit(d[keep, ])
  for (family in c("npqreg", "npconmode")) {
    fun <- get(family)
    f <- if (family == "npqreg") y ~ x else cy ~ x
    h <- if (family == "npqreg") c(.7, .6) else c(.2, .6)
    actual <- do.call(fun, list(formula = f, data = d, bws = h, subset = keep,
      na.action = na.omit, newdata = clean[1:5, ]))
    native <- fun(bws = actual$bws, txdat = clean["x"],
      tydat = if (family == "npqreg") clean$y else clean$cy,
      exdat = clean[1:5, "x", drop = FALSE])
    field <- if (family == "npqreg") "quantile" else "conmode"
    expect_equal(actual[[field]], native[[field]], tolerance = 1e-12)
  }
})
