test_that("formula constructors and one-call fits construct a formula once", {
  old <- options(np.messages = FALSE); on.exit(options(old), add = TRUE)
  set.seed(1717)
  d <- data.frame(y = rnorm(24), x = runif(24), z = rnorm(24),
                  cy = factor(rep(1:2, 12)))
  cases <- list(
    npregbw = list("y ~ x", .4), npudensbw = list("~ x", .4),
    npudistbw = list("~ x", .4), npcdensbw = list("y ~ x", c(.4,.4)),
    npcdistbw = list("y ~ x", c(.4,.4)),
    npplregbw = list("y ~ x | z", matrix(.4,2,1)),
    npscoefbw = list("y ~ x | z", .4),
    npindexbw = list("y ~ x + z", c(1,.5,.4)),
    npksum = list("~ x", .4),
    npreg = list("y ~ x", .4), npudens = list("~ x", .4),
    npudist = list("~ x", .4), npcdens = list("y ~ x", c(.4,.4)),
    npcdist = list("y ~ x", c(.4,.4)),
    npplreg = list("y ~ x | z", matrix(.4,2,1)),
    npscoef = list("y ~ x | z", .4),
    npindex = list("y ~ x + z", c(1,.5,.4)),
    npqreg = list("y ~ x", c(.4,.4)), npconmode = list("cy ~ x", c(.2,.4)))
  count <- new.env(parent = emptyenv()); count$n <- 0L
  make.formula <- function() {
    count$n <- count$n + 1L
    if (count$n != 1L) stop("formula constructed again")
    runif(1L)
    as.formula(spec[[1L]], env = environment())
  }
  fields <- function(z) {
    if (inherits(z, "bandwidth") || inherits(z, "rbandwidth") ||
        inherits(z, "conbandwidth") || inherits(z, "condbandwidth") ||
        inherits(z, "sibandwidth") || inherits(z, "scbandwidth"))
      return(z$bw)
    if (inherits(z, "plbandwidth")) return(lapply(z$bw, `[[`, "bw"))
    if (inherits(z, "npkernelsum")) return(z$ksum)
    fitted(z)
  }
  for (family in names(cases)) {
    fun <- get(family, mode = "function")
    spec <- cases[[family]]
    controls <- list(data = quote(d), bws = spec[[2L]])
    if (family != "npksum") controls$bandwidth.compute <- FALSE
    if (family == "npindex") controls$se <- FALSE
    for (named in c(TRUE, FALSE)) {
      args <- c(if (named) list(formula = quote(make.formula())) else
        list(quote(make.formula())), controls)
      set.seed(1718); count$n <- 0L
      actual <- eval(as.call(c(list(quote(fun)), args)))
      actual.rng <- .Random.seed
      expect_identical(count$n, 1L, info = paste(family, named))
      set.seed(1718); count$n <- 0L
      resolved <- make.formula()
      args[[1L]] <- resolved
      expected <- eval(as.call(c(list(quote(fun)), args)))
      expect_identical(actual.rng, .Random.seed, info = paste(family, named))
      expect_identical(fields(actual), fields(expected), info = paste(family, named))
    }
  }
})

test_that("formula handoff preserves environments, provenance and lazy subsets", {
  old <- options(np.messages = FALSE); on.exit(options(old), add = TRUE)
  d <- data.frame(y = sin(1:24), x = seq(.1,.9,length.out=24))
  count <- new.env(parent = emptyenv()); count$n <- 0L
  formula.env <- list2env(list(shift = .15), parent = baseenv())
  formula.env$shifted <- eval(quote(function(x) x + shift), formula.env)
  make.formula <- function() {
    count$n <- count$n + 1L
    if (count$n > 1L) stop("formula replay")
    as.formula("y ~ shifted(x)", env = formula.env)
  }
  wrapper <- function(...) npreg(...)
  actual <- wrapper(formula = make.formula(), data = d, subset = x > .2,
                    bws = .4, bandwidth.compute = FALSE)
  expect_identical(count$n, 1L)
  expect_identical(environment(actual$bws$formula), formula.env)
  expect_identical(actual$bws$call$formula, quote(make.formula()))
  expected <- npreg(txdat = data.frame(x = d$x[d$x > .2] + .15),
                    tydat = d$y[d$x > .2], bws = .4, bandwidth.compute = FALSE)
  expect_identical(fitted(actual), fitted(expected))
  expect_identical(fitted(npreg(actual$bws)), fitted(actual))
  expect_equal(predict(actual, newdata = d[1:4, ]),
               predict(expected, exdat = data.frame(x = d$x[1:4] + .15)), tolerance = 0)
  expect_identical(count$n, 1L)
})

test_that("a formula construction error is not retried or substituted", {
  count <- new.env(parent = emptyenv()); count$n <- 0L
  failure <- structure(list(message = "formula factory sentinel", call = NULL),
                       class = c("np_factory_test_error", "error", "condition"))
  broken <- function() { count$n <- count$n + 1L; stop(failure) }
  for (family in c("npregbw", "npcdensbw", "npplregbw", "npksum", "npcdens")) {
    count$n <- 0L
    result <- tryCatch(get(family)(formula = broken()), error = identity)
    expect_identical(result, failure)
    expect_identical(count$n, 1L)
  }
})

test_that("formula factories preserve reordered and forwarded call shapes", {
  old <- options(np.messages = FALSE); on.exit(options(old), add = TRUE)
  d <- data.frame(y = sin(1:24), x = seq(.1,.9,length.out=24))
  count <- 0L
  make.formula <- function() {
    count <<- count + 1L
    if (count > 1L) stop("formula replay")
    y ~ x
  }
  forwarded <- function(...) npcdens(...)
  nested <- function(...) forwarded(...)
  oracle <- npcdens(y ~ x, data = d, bws = c(.4,.4), bandwidth.compute = FALSE)
  calls <- list(
    quote(npcdens(data = d, bws = c(.4,.4), formula = make.formula(), bandwidth.compute = FALSE)),
    quote(npcdens(bws = c(.4,.4), formula = make.formula(), data = d, bandwidth.compute = FALSE)),
    quote(forwarded(make.formula(), data = d, bws = c(.4,.4), bandwidth.compute = FALSE)),
    quote(nested(formula = make.formula(), data = d, bws = c(.4,.4), bandwidth.compute = FALSE)))
  for (expr in calls) {
    count <- 0L
    actual <- eval(expr)
    expect_identical(count, 1L)
    expect_identical(fitted(actual), fitted(oracle))
    expect_identical(fitted(npcdens(actual$bws)), fitted(oracle))
    expect_identical(count, 1L)
  }
  # These inherited unsupported call shapes must not be reinterpreted as a
  # different formula/data route by the value handoff. They are not positive
  # coverage claims; separate named formula= is the supported reordered form.
  count <- 0L
  expect_error(npcdens(bws = c(.4,.4), make.formula(), data = d,
                       bandwidth.compute = FALSE), "xdat must be a data frame")
  expect_identical(count, 1L)
  count <- 0L
  expect_error(npqreg(bws = make.formula(), txdat = d["x"], tydat = d$y,
                      bwmethod = "normal-reference", tau = .5),
               "unused argument in npcdistbw: '.np.formula.state'", fixed = TRUE)
  expect_identical(count, 1L)
})
