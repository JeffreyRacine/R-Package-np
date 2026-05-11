test_that("npconmode proper helper enforces binary complement probabilities", {
  helper <- getFromNamespace(".npConmodeProperProbabilities", "np")
  raw <- matrix(c(-0.2, 1.2,
                  0.4, 0.8,
                  0.3, 0.3), ncol = 2L, byrow = TRUE)

  out <- helper(raw, levels = c("0", "1"), proper = TRUE)

  expect_true(isTRUE(out$proper.requested))
  expect_true(isTRUE(out$proper.applied))
  expect_equal(rowSums(out$probabilities), rep(1, nrow(raw)), tolerance = 1e-12)
  expect_true(all(out$probabilities >= -1e-12))
  expect_equal(out$probabilities[, 2L], 1 - out$probabilities[, 1L], tolerance = 1e-12)
  expect_identical(out$proper.info$reason, "projected")
})

test_that("npconmode proper controls fail fast on invalid control objects", {
  helper <- getFromNamespace(".npConmodeProperProbabilities", "np")
  raw <- matrix(c(0.2, 0.8, 0.4, 0.6), ncol = 2L, byrow = TRUE)

  expect_error(helper(raw, levels = c("0", "1"), proper.control = 1),
               "'proper.control' must be a list", fixed = TRUE)
  expect_error(helper(raw, levels = c("0", "1"), proper.control = list(foo = 1)),
               "unused argument")
  expect_error(helper(raw, levels = c("0", "1"), proper.control = list(tol = -1)),
               "'proper.control$tol' must be a non-negative scalar", fixed = TRUE)
})

test_that("npconmode bandwidth route reports all-NA training data clearly", {
  set.seed(20260511)
  d <- data.frame(x = runif(30), y = factor(rbinom(30, 1L, .5), levels = 0:1))
  bw <- npcdensbw(y ~ x, data = d, nmulti = 1L)
  badx <- data.frame(x = rep(NA_real_, 10))
  bady <- data.frame(y = factor(rep(NA, 10), levels = 0:1))

  expect_error(npconmode(bws = bw, txdat = badx, tydat = bady),
               "Training data has no rows without NAs", fixed = TRUE)
})

test_that("npconmode proper defaults follow the canonical regression type", {
  effective <- getFromNamespace(".npConmodeEffectiveProper", "np")

  expect_false(effective(list(regtype = "lc"), NULL))
  expect_true(effective(list(regtype = "ll"), NULL))
  expect_true(effective(list(regtype = "lp"), NULL))
  expect_false(effective(list(regtype = "lp", regtype.engine = "lc"), NULL))
  expect_true(effective(list(regtype = "lc", regtype.engine = "lp"), NULL))
  expect_false(effective(list(regtype = "lp"), FALSE))
  expect_true(effective(list(regtype = "lc"), TRUE))
})

test_that("npconmode formula route forwards NOMAD shortcut controls", {
  skip_if_not_installed("crs")

  set.seed(20260511)
  n <- 30L
  d <- data.frame(
    x = seq(-1, 1, length.out = n)
  )
  d$y <- factor(rbinom(n, 1L, plogis(1.5 * d$x)))

  capture.output(
    fit <- npconmode(
      y ~ x,
      data = d,
      nomad = TRUE,
      nmulti = 1L,
      nomad.nmulti = 1L,
      degree.max = 1L
    )
  )

  expect_s3_class(fit, "conmode")
  expect_true(isTRUE(fit$bws$nomad.shortcut$enabled))
  expect_identical(as.character(fit$bws$regtype.engine), "lp")
})

test_that("npconmode NOMAD fit metadata propagates through evaluation and bandwidth plots", {
  skip_if_not_installed("crs")

  count_namespace_calls <- function(fun, expr) {
    ns <- asNamespace("np")
    ctr <- new.env(parent = emptyenv())
    ctr$n <- 0L
    trace(
      what = fun,
      where = ns,
      tracer = bquote(.(ctr)$n <- .(ctr)$n + 1L),
      print = FALSE
    )
    on.exit(try(untrace(fun, where = ns), silent = TRUE), add = TRUE)
    force(expr)
    ctr$n
  }

  set.seed(20260511)
  n <- 36L
  d <- data.frame(x = seq(-1, 1, length.out = n))
  d$y <- factor(rbinom(n, 1L, plogis(1.5 * d$x)))
  nd <- data.frame(x = c(-0.75, -0.1, 0.35, 0.8))

  capture.output(
    fit <- npconmode(
      y ~ x,
      data = d,
      newdata = nd,
      nomad = TRUE,
      nmulti = 1L,
      nomad.nmulti = 1L,
      degree.max = 1L,
      probabilities = TRUE
    )
  )

  expect_s3_class(fit, "conmode")
  expect_true(isTRUE(fit$bws$nomad.shortcut$enabled))
  expect_identical(as.character(fit$bws$regtype.engine), "lp")
  expect_length(fit$conmode, nrow(nd))
  expect_length(fitted(fit), nrow(nd))
  expect_equal(nrow(fit$probabilities), nrow(nd))
  expect_true(all(fit$probabilities >= -1e-12))
  expect_equal(rowSums(fit$probabilities), rep(1, nrow(nd)), tolerance = 1e-12)

  bw.calls <- count_namespace_calls("npcdensbw", {
    out <- plot(
      fit$bws,
      xdat = d["x"],
      ydat = d["y"],
      output = "data",
      perspective = FALSE,
      neval = 6L
    )
    expect_type(out, "list")
  })
  expect_identical(bw.calls, 0L)
})

test_that("npconmode binary class-probability gradients are stored and plotted", {
  set.seed(20260511)
  n <- 70L
  d <- data.frame(x = seq(-1, 1, length.out = n))
  d$y <- factor(rbinom(n, 1L, plogis(1.25 * d$x)), levels = 0:1)

  capture.output(
    fit <- npconmode(
      y ~ x,
      data = d,
      regtype = "ll",
      bwmethod = "cv.ls",
      nmulti = 1L,
      probabilities = TRUE,
      gradients = TRUE
    )
  )

  expect_s3_class(fit, "conmode")
  expect_true(isTRUE(fit$gradients))
  expect_equal(rowSums(fit$probabilities), rep(1, n), tolerance = 1e-10)
  expect_equal(as.character(fit$probability.gradient.level), "0")
  expect_equal(predict(fit, type = "class"), fit$conmode)
  expect_equal(predict(fit, type = "prob"), fit$probabilities)

  gout <- gradients(fit)
  expect_equal(dim(gout), c(n, 1L))

  direct <- npcdens(
    bws = fit$bws,
    txdat = d["x"],
    tydat = d["y"],
    exdat = d["x"],
    eydat = factor(rep("0", n), levels = levels(d$y)),
    gradients = TRUE
  )
  expect_equal(as.vector(gout), as.vector(direct$congrad), tolerance = 1e-10)
  expect_error(gradients(fit, errors = TRUE), "standard errors")
  expect_error(gradients(fit, level = "1"), "stored class-probability gradients")

  pdata <- plot(fit, output = "data")
  expect_named(pdata, "x")
  expect_equal(nrow(pdata$x), n)

  grdata <- plot(fit, gradients = TRUE, output = "data")
  expect_named(grdata, "x")
  expect_true("effect" %in% names(grdata$x))
  expect_error(plot(fit, gradients = TRUE, level = "1", output = "data"),
               "stored class-probability gradients")

  pdf(tempfile(fileext = ".pdf"))
  expect_silent(plot(fit))
  expect_silent(plot(fit, gradients = TRUE))
  expect_error(plot(fit, errors = "bootstrap"), "unused plot argument")
  grDevices::dev.off()
})

test_that("npconmode class-probability gradients are level-specific for multinomial responses", {
  set.seed(20260511)
  n <- 45L
  d <- data.frame(x = seq(-1, 1, length.out = n))
  d$y <- factor(sample(letters[1:3], n, replace = TRUE))

  capture.output(
    fit <- npconmode(
      y ~ x,
      data = d,
      regtype = "ll",
      bwmethod = "cv.ls",
      nmulti = 1L,
      probabilities = TRUE,
      level = "b",
      gradients = TRUE
    )
  )

  expect_s3_class(fit, "conmode")
  expect_equal(as.character(fit$probability.gradient.level), "b")
  expect_equal(dim(gradients(fit)), c(n, 1L))
  expect_equal(rowSums(fit$probabilities), rep(1, n), tolerance = 1e-10)
  expect_error(gradients(fit, level = "c"), "stored class-probability gradients")
  grdata <- plot(fit, gradients = TRUE, output = "data")
  expect_equal(unique(grdata$x$level), "b")
})
