test_that("unconditional replay resolves data in its retained call owner", {
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  fixture <- data.frame(x = seq(-1, 1, length.out = 16L), z = sin(1:16))
  f <- ~ x + z
  for (family in c("npudens", "npudist")) {
    bw.fun <- get(paste0(family, "bw"))
    fit.fun <- get(family)
    build <- function() {
      private.data <- fixture
      bw.fun(f, data = private.data, bws = c(.4, .3), bandwidth.compute = FALSE)
    }
    bw <- unserialize(serialize(build(), NULL))
    reference <- fit.fun(bw, tdat = fixture)
    expect_identical(environment(fit.fun(bw)$bws$call), environment(bw$call))
    expect_identical(fitted(fit.fun(bw)), fitted(reference))
    expect_identical(fitted(fit.fun(bw, data = NULL)), fitted(reference))
    expect_identical(predict(fit.fun(bw), newdata = fixture[1:3, ]),
                     predict(reference, newdata = fixture[1:3, ]))
    changed <- transform(fixture, x = x * 2)
    assign("private.data", changed, envir = environment(bw$call))
    expect_identical(fitted(fit.fun(bw)), fitted(fit.fun(bw, data = changed)))
    expect_identical(fitted(fit.fun(bw, data = fixture)), fitted(reference))
    rm("private.data", envir = environment(bw$call))
    expect_error(fit.fun(bw), "private.data.*not found")
    expect_identical(fitted(fit.fun(bw, data = fixture)), fitted(reference))
  }
})

test_that("unconditional formula fits share one prepared training sample", {
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  set.seed(847L)
  d <- data.frame(x = rnorm(30), z = rnorm(30))
  for (family in c("npudens", "npudist")) {
    fit.fun <- get(family)
    count <- 0L
    counted <- function(x) { count <<- count + 1L; x }
    f <- ~ counted(x) + z
    fit <- fit.fun(f, data = d, bws = c(.5, .6))
    expect_identical(count, 1L)
    count <- 0L
    named <- fit.fun(formula = f, data = d, bws = c(.5, .6))
    expect_identical(count, 1L)
    expect_identical(fitted(fit), fitted(named))
    count <- 0L
    refit <- fit.fun(bws = fit$bws)
    expect_identical(count, 1L)
    expect_identical(fitted(fit), fitted(refit))
    count <- 0L
    auto <- fit.fun(f, data = d, bwmethod = "normal-reference", se = TRUE)
    expect_identical(count, 1L)
    expect_true(isTRUE(auto[["se", exact = TRUE]]))
    expect_false(".np.formula.state" %in% names(auto$bws$call))
    expect_false(".np.formula.state" %in% names(auto$bws$call$...))
    count <- 0L
    named.auto <- fit.fun(formula = f, data = d, bwmethod = "normal-reference", se = TRUE)
    expect_identical(count, 1L)
    expect_identical(fitted(named.auto), fitted(auto))
    expect_error(fit.fun(formula = ~ .np_missing_function(x), data = d,
                        bws = .5), "could not find function")
    recovered <- fit.fun(formula = f, data = d, bws = c(.5, .6))
    expect_identical(fitted(recovered), fitted(fit))
    count <- 0L
    new <- fit.fun(bws = fit$bws, newdata = d[1:4, , drop = FALSE], se = TRUE)
    expect_identical(count, 2L)
    control <- fit.fun(bws = fit$bws, tdat = d, edat = d[1:4, , drop = FALSE], se = TRUE)
    expect_equal(fitted(new), fitted(control), tolerance = 1e-14)
    expect_equal(se(new), se(control), tolerance = 1e-14)
    jittered <- function(x) x + runif(length(x), -.1, .1)
    set.seed(848L)
    actual <- fit.fun(~ jittered(x) + z, data = d, bws = c(.5, .6))
    rng <- .Random.seed
    set.seed(848L)
    native <- d
    native$x <- jittered(d$x)
    expected <- fit.fun(tdat = native, bws = c(.5, .6))
    expect_identical(.Random.seed, rng)
    expect_equal(fitted(actual), fitted(expected), tolerance = 1e-14)
  }
})
