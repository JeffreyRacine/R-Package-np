test_that("legacy conditional alignment is unwrapped without replacing user prediction metadata", {
  unwrap <- get(".np_formula_unwrap_prediction", envir = environment(npregbw))
  tt <- terms(~ y + x)
  variables <- attr(tt, "variables")
  attr(tt, "predvars") <- quote(as.data.frame(ts.intersect(y, x)))
  expect_identical(unwrap(tt)$prediction, variables)
  expect_true(unwrap(tt)$makepredictcall)
  for (user in list(quote(list(y, scale(x, center = 2, scale = 3))),
                    quote(as.data.frame(custom(y, x))),
                    quote(as.data.frame(ts.intersect(x, y))))) {
    attr(tt, "predvars") <- user
    expect_identical(unwrap(tt)$prediction, user)
    expect_false(unwrap(tt)$makepredictcall)
  }
  tt <- terms(~ u + y + x)
  variables <- attr(tt, "variables")
  attr(tt, "predvars") <- substitute(cbind(as.data.frame(ts.intersect(y, x)), u,
    check.rows = TRUE)[, INDEX], list(INDEX = c(3L, 1L, 2L)))
  expect_identical(unwrap(tt)$prediction, variables)
  expect_true(unwrap(tt)$makepredictcall)
  user <- substitute(cbind(as.data.frame(ts.intersect(y, x)), u,
    check.rows = TRUE)[, INDEX], list(INDEX = c(1L, 2L, 3L)))
  attr(tt, "predvars") <- user
  expect_identical(unwrap(tt)$prediction, user)
  expect_false(unwrap(tt)$makepredictcall)
  tt <- terms(~ stop("never evaluate while identifying the scaffold") + x)
  original <- attr(tt, "variables")
  attr(tt, "predvars") <- as.call(list(quote(as.data.frame),
    as.call(c(list(quote(ts.intersect)), as.list(original)[-1L]))))
  expect_identical(unwrap(tt)$prediction, original)
})

test_that("conditional fitting reuses one prepared training frame", {
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  set.seed(865L)
  d <- data.frame(y = rnorm(48), x = rnorm(48), z = factor(rep(1:2, 24)))
  counter <- new.env(parent = emptyenv())
  counter$n <- 0L
  counted <- function(x) { counter$n <- counter$n + 1L; x }
  for (family in c("npcdens", "npcdist")) {
    bw.fun <- get(paste0(family, "bw"))
    fit.fun <- get(family)
    f <- y ~ counted(x) + z
    counter$n <- 0L
    bw <- bw.fun(f, data = d, bws = c(.65, .55, .2), bandwidth.compute = FALSE)
    expect_identical(counter$n, 1L)
    counter$n <- 0L
    fit <- fit.fun(f, data = d, bws = c(.65, .55, .2))
    expect_identical(counter$n, 1L)
    counter$n <- 0L
    named <- fit.fun(formula = f, data = d, bws = c(.65, .55, .2))
    expect_identical(counter$n, 1L)
    expect_identical(fitted(named), fitted(fit))
    counter$n <- 0L
    refit <- fit.fun(bws = fit$bws)
    expect_identical(counter$n, 0L)
    expect_identical(fitted(refit), fitted(fit))
    counter$n <- 0L
    auto <- fit.fun(formula = f, data = d, bandwidth.compute = FALSE)
    expect_identical(counter$n, 1L)
    expect_false(".np.formula.state" %in% names(auto$bws$call$...))
    expect_false(grepl("getFromNamespace",
      paste(deparse(attr(terms(bw), "predvars")), collapse = ""), fixed = TRUE))
    counter$n <- 0L
    evaluated <- fit.fun(bws = bw, newdata = d[1:5, ], se = TRUE, gradients = TRUE)
    expect_identical(counter$n, 1L)
    native <- fit.fun(bws = bw, txdat = d[c("x", "z")], tydat = d["y"],
      exdat = d[1:5, c("x", "z")], eydat = d[1:5, "y", drop = FALSE],
      se = TRUE, gradients = TRUE)
    expect_equal(fitted(evaluated), fitted(native), tolerance = 1e-12)
    expect_equal(se(evaluated), se(native), tolerance = 1e-12)
    expect_equal(unname(gradients(evaluated)), unname(gradients(native)), tolerance = 1e-12)
    jittered <- function(x) x + runif(length(x), -.1, .1)
    set.seed(868L)
    actual <- fit.fun(y ~ jittered(x) + z, data = d, bws = c(.65, .55, .2))
    rng <- .Random.seed
    set.seed(868L)
    dx <- d[c("x", "z")]; dx$x <- jittered(dx$x)
    expected <- fit.fun(txdat = dx, tydat = d["y"], bws = c(.65, .55, .2))
    expect_identical(.Random.seed, rng)
    expect_equal(fitted(actual), fitted(expected), tolerance = 1e-12)
  }
})

test_that("conditional formulas retain aligned roles and trained prediction terms", {
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  set.seed(869L)
  y <- ts(rnorm(40), frequency = 4)
  z <- factor(rep(c("a", "b", "c"), length.out = 39))
  x <- data.frame(a = as.numeric(y)[1:39], z = z)
  response <- as.numeric(y)[2:40]
  for (family in c("npcdens", "npcdist")) {
    bw.fun <- get(paste0(family, "bw")); fit.fun <- get(family)
    fit <- fit.fun(y ~ lag(y, -1) + z, bws = c(.7, .8, .2), se = TRUE)
    native <- fit.fun(bws = fit$bws, txdat = x, tydat = response, se = TRUE)
    expect_length(fitted(fit), 39L)
    expect_equal(fitted(fit), fitted(native), tolerance = 1e-12)
    expect_equal(se(fit), se(native), tolerance = 1e-12)
    d <- data.frame(y = response, x = x$a, z = z)
    bw <- bw.fun(y ~ log(abs(x) + 1) + z, data = d, bws = c(.7, .8, .2),
      bandwidth.compute = FALSE)
    new <- transform(d[1:5, ], x = x + 2)
    evaluated <- fit.fun(bws = bw, newdata = new)
    sx <- log(abs(d$x) + 1)
    ex <- log(abs(new$x) + 1)
    expected <- fit.fun(bws = bw, txdat = data.frame(x = as.numeric(sx), z = z),
      tydat = response, exdat = data.frame(x = ex, z = new$z), eydat = new["y"])
    expect_equal(fitted(evaluated), fitted(expected), tolerance = 1e-12)
    # Scalar-matrix fitting remains outside this family's existing type
    # contract. The common frame owner must still retain prediction metadata.
    sbw <- bw.fun(y ~ scale(x) + z, data = d, bws = c(.7, .8, .2),
      bandwidth.compute = FALSE)
    frame <- get(".np_formula_model_frame", envir = environment(bw.fun))
    scaled <- scale(d$x)
    prepared <- frame(terms(sbw), data = new)
    expect_equal(as.numeric(prepared[["scale(x)"]]),
      (new$x - attr(scaled, "scaled:center")) / attr(scaled, "scaled:scale"),
      tolerance = 1e-14)
    expect_error(fit.fun(sbw), "supplied bandwidths do not match 'txdat' in type")
    build <- function() {
      private.data <- d
      bw.fun(y ~ x + z, data = private.data, bws = c(.7, .8, .2),
        bandwidth.compute = FALSE)
    }
    saved <- unserialize(serialize(build(), NULL))
    expect_equal(fitted(fit.fun(saved)),
      fitted(fit.fun(saved, txdat = d[c("x", "z")], tydat = d["y"])),
      tolerance = 1e-12)
  }
})
