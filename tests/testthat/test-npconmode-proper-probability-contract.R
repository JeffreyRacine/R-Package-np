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

test_that("npconmode rejects non-categorical responses before modal fitting", {
  set.seed(20260513)
  d <- data.frame(x = runif(30), y = rnorm(30))

  expect_error(
    npconmode(y ~ x, data = d, bws = c(0.4, 0.4),
              bandwidth.compute = FALSE),
    "'tydat' must consist of one (1) discrete variable",
    fixed = TRUE
  )
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

test_that("npconmode formula route treats response-free newdata strictly", {
  nested_call <- function() {
    set.seed(20260511)
    d <- data.frame(x = sort(runif(50)))
    d$y <- factor(rbinom(50, 1L, plogis(1.2 * d$x)), levels = 0:1)
    nd <- data.frame(x = c(.1, .4, .8))

    npconmode(
      y ~ x,
      data = d,
      newdata = nd,
      regtype = "ll",
      bwmethod = "cv.ls",
      nmulti = 1L,
      probabilities = TRUE,
      gradients = TRUE,
      level = "1"
    )
  }

  fit <- nested_call()
  expect_s3_class(fit, "conmode")
  expect_equal(length(fitted(fit)), 3L)
  expect_equal(nrow(predict(fit, type = "prob")), 3L)
  expect_equal(NROW(gradients(fit)), 3L)
})

test_that("npconmode formula route keeps native exdat out of bandwidth selection", {
  set.seed(20260514)
  d <- data.frame(x = sort(runif(40)))
  d$y <- factor(rbinom(40, 1L, plogis(0.5 + d$x)), levels = 0:1)
  nd <- data.frame(x = c(0.2, 0.5, 0.8))

  fit <- npconmode(
    y ~ x,
    data = d,
    exdat = nd,
    probabilities = TRUE,
    nmulti = 1
  )

  expect_s3_class(fit, "conmode")
  expect_equal(length(fitted(fit)), nrow(nd))
  expect_equal(nrow(predict(fit, type = "prob")), nrow(nd))
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
  expect_equal(dim(fit$probability.errors), dim(fit$probabilities))
  expect_equal(colnames(fit$probability.errors), colnames(fit$probabilities))
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
  expect_equal(plot(fit, output = "both"), pdata)

  perr <- plot(fit, output = "data", errors = "asymptotic",
               band = "all", level = "0")
  expect_true(all(c("stderr", "lower", "upper",
                    "pointwise.lower", "bonferroni.upper") %in%
                    names(perr$x)))
  expect_equal(perr$x$stderr, fit$probability.errors[, "0"],
               tolerance = 1e-12, ignore_attr = TRUE)

  grdata <- plot(fit, gradients = TRUE, output = "data")
  expect_named(grdata, "x")
  expect_true("effect" %in% names(grdata$x))
  expect_error(plot(fit, gradients = TRUE, level = "1", output = "data"),
               "stored class-probability gradients")

  grid.data <- plot(fit, view = "fixed", neval = 9L, output = "data")
  expect_named(grid.data, "x")
  expect_equal(nrow(grid.data$x), 9L)
  expect_equal(unique(grid.data$x$view), "fixed")
  grid.prob <- predict(fit,
                       newdata = data.frame(x = grid.data$x$x),
                       type = "prob")[, "0"]
  expect_equal(grid.data$x$probability, grid.prob, tolerance = 1e-5)

  grid.err <- plot(fit, view = "fixed", neval = 9L, output = "data",
                   errors = "asymptotic", band = "pointwise", level = "0")
  expect_true(all(c("stderr", "lower", "upper") %in% names(grid.err$x)))
  expect_true(any(is.finite(grid.err$x$stderr)))
  pred.se <- predict(fit,
                     newdata = data.frame(x = grid.data$x$x),
                     type = "prob",
                     se.fit = TRUE)
  expect_named(pred.se, c("fit", "se.fit"))
  expect_equal(pred.se$fit[, "0"], grid.prob, tolerance = 1e-5)
  expect_true(any(is.finite(pred.se$se.fit[, "0"])))
  expect_error(predict(fit, type = "class", se.fit = TRUE),
               "type=\"prob\"", fixed = TRUE)

  set.seed(20260514)
  grid.boot <- plot(fit, view = "fixed", neval = 7L, output = "data",
                    errors = "bootstrap", bootstrap = "inid", B = 39L,
                    band = "pointwise", level = "0")
  expect_true(all(c("stderr", "center", "lower", "upper") %in%
                    names(grid.boot$x)))
  expect_true(any(is.finite(grid.boot$x$stderr)))
  expect_true(all(is.finite(grid.boot$x$lower)))
  expect_true(all(is.finite(grid.boot$x$upper)))
  expect_error(plot(fit, view = "fixed", neval = 7L, output = "data",
                    errors = "bootstrap", bootstrap = "wild"),
               "wild bootstrap is not available")

  grid.grad <- plot(fit, view = "fixed", neval = 9L,
                    gradients = TRUE, output = "data")
  expect_equal(nrow(grid.grad$x), 9L)
  expect_true("effect" %in% names(grid.grad$x))

  pdf(tempfile(fileext = ".pdf"))
  expect_silent(plot(fit))
  expect_silent(plot(fit, view = "fixed", neval = 9L))
  expect_silent(plot(fit, col = "red", lwd = 2, lty = 2, type = "l"))
  expect_silent(plot(fit, gradients = TRUE))
  expect_silent(plot(fit, gradients = TRUE, view = "fixed", neval = 9L))
  expect_silent(plot(fit, gradients = TRUE, col = "blue", lwd = 2, lty = 3,
                     type = "l"))
  expect_silent(plot(fit, errors = "bootstrap", B = 39L))
  expect_error(plot(fit, gradients = TRUE, errors = "asymptotic"),
               "not probability gradients")
  expect_error(plot(fit, renderer = "rgl"),
               "renderer='rgl' is supported only")
  grDevices::dev.off()
})

test_that("npconmode fixed surface payload is object-fed and predict-aligned", {
  set.seed(20260514)
  n <- 60L
  d <- data.frame(
    x1 = runif(n, -1, 1),
    x2 = runif(n, -1, 1),
    z = factor(sample(c("a", "b"), n, replace = TRUE))
  )
  eta <- 0.5 + d$x1 - 0.75 * d$x2 + 0.25 * (d$z == "b")
  d$y <- factor(rbinom(n, 1L, plogis(eta)), levels = 0:1)

  capture.output(
    fit <- npconmode(
      y ~ x1 + x2 + z,
      data = d,
      regtype = "ll",
      bwmethod = "cv.ls",
      nmulti = 1L,
      probabilities = TRUE
    )
  )

  surf <- plot(fit, view = "fixed", perspective = TRUE,
               neval = 7L, output = "data")
  expect_s3_class(surf, "np_conmode_surface_data")
  expect_named(surf, c("surface", "grid", "held", "level", "levels", "proper"))
  expect_equal(surf$grid$variables, c("x1", "x2"))
  expect_equal(surf$grid$dims, c(7L, 7L))
  expect_equal(nrow(surf$surface), 49L)
  expect_true("z" %in% names(surf$held))
  expect_equal(unique(surf$surface$view), "fixed")
  expect_true(all(is.finite(surf$surface$probability)))
  expect_true(all(surf$surface$probability >= -1e-10))
  expect_true(all(surf$surface$probability <= 1 + 1e-10))

  surf.err <- plot(fit, view = "fixed", perspective = TRUE,
                   neval = 5L, output = "data",
                   errors = "asymptotic", level = surf$level)
  expect_true(all(c("stderr", "lower", "upper") %in% names(surf.err$surface)))
  expect_true(any(is.finite(surf.err$surface$stderr)))
  expect_error(plot(fit, view = "fixed", perspective = TRUE,
                    neval = 5L, errors = "asymptotic"),
               "output=\"data\" only", fixed = TRUE)
  expect_error(plot(fit, view = "fixed", perspective = TRUE,
                    neval = 5L, errors = "bootstrap", output = "data",
                    B = 39L),
               "surface bootstrap intervals")

  nd <- data.frame(
    x1 = surf$surface$x1,
    x2 = surf$surface$x2,
    z = surf$held$z[rep(1L, nrow(surf$surface))]
  )
  pred <- predict(fit, newdata = nd, type = "prob")[, surf$level]
  expect_equal(surf$surface$probability, pred, tolerance = 1e-5)

  surf.alias <- plot(fit, view = "fixed", persp = TRUE,
                     neval = 5L, output = "data",
                     plot.vars = c("x2", "x1"))
  expect_equal(surf.alias$grid$variables, c("x2", "x1"))
  expect_equal(surf.alias$grid$dims, c(5L, 5L))

  expect_silent({
    grDevices::pdf(tempfile(fileext = ".pdf"))
    on.exit(grDevices::dev.off(), add = TRUE)
    plot(fit, perspective = TRUE, neval = 5L)
    rendered <- plot(fit, perspective = TRUE, neval = 5L,
                     output = "plot-data")
  })
  data.only <- plot(fit, perspective = TRUE, neval = 5L,
                    output = "data")
  expect_equal(rendered, data.only)
  expect_equal(
    plot(fit, perspective = TRUE, renderer = "rgl", neval = 5L,
         output = "data"),
    data.only
  )
  expect_error(plot(fit, renderer = "rgl", output = "data"),
               "renderer='rgl' is supported only")
  if (suppressWarnings(requireNamespace("rgl", quietly = TRUE))) {
    rgl.rendered <- suppressWarnings(
      plot(fit, perspective = TRUE, renderer = "rgl", neval = 5L,
           output = "plot-data")
    )
    expect_equal(rgl.rendered, data.only)
  }
  expect_error(plot(fit, plot.vars = c("x1", "x2"), output = "data"),
               "'plot.vars' is only supported")
  expect_error(plot(fit, perspective = TRUE, persp = FALSE, output = "data"),
               "conflicting plot.conmode arguments")
  expect_error(plot(fit, perspective = TRUE, gradients = TRUE, output = "data"),
               "not gradients")
  expect_error(plot(fit, perspective = TRUE, plot.vars = c("x1", "z"),
                    output = "data"),
               "currently require continuous variables")

  fit.one <- npconmode(
    y ~ x1 + z,
    data = d,
    bws = c(0.5, 0.5, 0.5),
    bandwidth.compute = FALSE,
    probabilities = TRUE
  )
  expect_error(plot(fit.one, perspective = TRUE, output = "data"),
               "require exactly two continuous conditioning variables")
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
  expect_true(all(fit$probabilities >= -1e-12))
  expect_equal(rowSums(fit$probabilities), rep(1, n), tolerance = 1e-10)
  expect_error(
    npconmode(
      y ~ x,
      data = d,
      regtype = "ll",
      bwmethod = "cv.ls",
      nmulti = 1L,
      probabilities = TRUE,
      level = "z",
      gradients = TRUE
    ),
    "'level' must identify one response level",
    fixed = TRUE
  )
  expect_error(gradients(fit, level = "c"), "stored class-probability gradients")
  grdata <- plot(fit, gradients = TRUE, output = "data")
  expect_equal(unique(grdata$x$level), "b")
})
