test_that("fixed scalar npcdens CVLS uses the canonical constant-X limit", {
  old <- options(
    np.messages = FALSE,
    np.tree = FALSE,
    np.categorical.compress = FALSE,
    np.largeh = FALSE,
    np.objective.cache = FALSE,
    np.largelambda = TRUE
  )
  on.exit(options(old), add = TRUE)
  eval.only <- getFromNamespace(".npcdensbw_eval_only", "npRmpi")

  set.seed(202608094L)
  n <- 120L
  y <- data.frame(y = rnorm(n))
  xu <- data.frame(z = factor(sample(letters[1:3], n, replace = TRUE)))
  xo <- data.frame(z = ordered(sample(1:4, n, replace = TRUE)))

  evaluate <- function(x, lambda, enabled) {
    bw <- npcdensbw(
      xdat = x,
      ydat = y,
      bws = c(0.3, lambda),
      bandwidth.compute = FALSE,
      bwmethod = "cv.ls",
      regtype = "lc"
    )
    options(np.largelambda = enabled)
    eval.only(x, y, bw)
  }

  unordered.on <- evaluate(xu, 2 / 3, TRUE)
  unordered.off <- evaluate(xu, 2 / 3, FALSE)
  unordered.on.again <- evaluate(xu, 2 / 3, TRUE)
  ordered.on <- evaluate(xo, 1, TRUE)
  ordered.off <- evaluate(xo, 1, FALSE)

  expect_equal(unordered.on$objective, unordered.off$objective,
               tolerance = 1e-10)
  expect_equal(unordered.on$objective, unordered.on.again$objective,
               tolerance = 1e-10)
  expect_equal(ordered.on$objective, unordered.on$objective,
               tolerance = 1e-10)
  expect_true(is.finite(ordered.off$objective))
  expect_equal(unordered.on$num.feval.fast, 1)
  expect_equal(unordered.off$num.feval.fast, 0)
  expect_equal(unordered.on.again$num.feval.fast, 1)
  expect_equal(ordered.on$num.feval.fast, 1)
  expect_equal(ordered.off$num.feval.fast, 0)
})

test_that("nonfixed conditional Y bandwidths do not enter the fixed X limit", {
  old <- options(
    np.messages = FALSE,
    np.tree = FALSE,
    np.categorical.compress = FALSE,
    np.largeh = FALSE,
    np.objective.cache = FALSE,
    np.largelambda = TRUE
  )
  on.exit(options(old), add = TRUE)
  eval.only <- getFromNamespace(".npcdensbw_eval_only", "npRmpi")

  set.seed(202608095L)
  n <- 100L
  x <- data.frame(z = factor(sample(letters[1:3], n, replace = TRUE)))
  y <- data.frame(y = rnorm(n))

  for (type in c("generalized_nn", "adaptive_nn")) {
    bw <- npcdensbw(
      xdat = x,
      ydat = y,
      bws = c(20, 2 / 3),
      bandwidth.compute = FALSE,
      bwmethod = "cv.ls",
      bwtype = type,
      regtype = "lc"
    )
    options(np.largelambda = TRUE)
    enabled <- eval.only(x, y, bw)
    options(np.largelambda = FALSE)
    disabled <- eval.only(x, y, bw)

    expect_equal(enabled$objective, disabled$objective, tolerance = 1e-10,
                 info = type)
    expect_equal(enabled$num.feval.fast, 0, info = type)
    expect_equal(disabled$num.feval.fast, 0, info = type)
  }
})

test_that("constant X remains side-specific when Y is categorical or bounded", {
  old <- options(
    np.messages = FALSE,
    np.tree = FALSE,
    np.categorical.compress = FALSE,
    np.largeh = FALSE,
    np.objective.cache = FALSE,
    np.largelambda = TRUE
  )
  on.exit(options(old), add = TRUE)
  eval.only <- getFromNamespace(".npcdensbw_eval_only", "npRmpi")

  set.seed(202608096L)
  n <- 120L
  x <- data.frame(x = factor(sample(letters[1:3], n, replace = TRUE)))

  y.cat <- data.frame(y = factor(sample(LETTERS[1:4], n, replace = TRUE)))
  bw.cat <- npcdensbw(
    xdat = x,
    ydat = y.cat,
    bws = c(0.2, 2 / 3),
    bandwidth.compute = FALSE,
    bwmethod = "cv.ls",
    regtype = "lc"
  )
  options(np.largelambda = TRUE)
  cat.enabled <- eval.only(x, y.cat, bw.cat)
  options(np.largelambda = FALSE)
  cat.disabled <- eval.only(x, y.cat, bw.cat)
  expect_equal(cat.enabled$objective, cat.disabled$objective,
               tolerance = 1e-10)
  expect_equal(cat.enabled$num.feval.fast, 1)
  expect_equal(cat.disabled$num.feval.fast, 0)

  y.beta <- data.frame(y = rbeta(n, 2, 3))
  bw.beta <- npcdensbw(
    xdat = x,
    ydat = y.beta,
    bws = c(0.2, 2 / 3),
    bandwidth.compute = FALSE,
    bwmethod = "cv.ls",
    regtype = "lc",
    cykertype = "beta",
    cykerbound = "range"
  )
  options(np.largelambda = TRUE)
  beta.enabled <- eval.only(x, y.beta, bw.beta)
  options(np.largelambda = FALSE)
  beta.disabled <- eval.only(x, y.beta, bw.beta)
  expect_equal(beta.enabled$objective, beta.disabled$objective,
               tolerance = 1e-10)
  expect_equal(beta.enabled$num.feval.fast, 0)
  expect_equal(beta.disabled$num.feval.fast, 0)
})
