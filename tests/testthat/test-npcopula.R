test_that("npcopula basic functionality works", {
  data("faithful")
  # npcopula can take a dbandwidth (from npudistbw) for copula 
  # or bandwidth (from npudensbw) for copula density
  
  bw <- npudistbw(dat=faithful, bws=c(0.5, 5), bandwidth.compute=FALSE)
  
  # Copula CDF
  cop <- npcopula(data=faithful, bws=bw)
  expect_s3_class(cop, "npcopula")
  expect_true(is.data.frame(as.data.frame(cop)))
  expect_true("copula" %in% names(cop))
  expect_identical(cop$target, "distribution")
  expect_identical(cop$evaluation, "sample")
  expect_true(all(cop$copula >= 0 & cop$copula <= 1))
  
  # Copula density
  bw_dens <- npudensbw(dat=faithful, bws=c(0.5, 5), bandwidth.compute=FALSE)
  cop_dens <- npcopula(data=faithful, bws=bw_dens)
  expect_s3_class(cop_dens, "npcopula")
  expect_true(is.data.frame(as.data.frame(cop_dens)))
  expect_true("copula" %in% names(cop_dens))
  expect_identical(cop_dens$target, "density")
  expect_identical(cop_dens$evaluation, "sample")
  expect_true(all(cop_dens$copula >= 0))
})

test_that("npcopula formula route builds a plot-ready grid by default", {
  set.seed(42)
  dat <- data.frame(x = rnorm(50), y = 0.5 * rnorm(50) + rnorm(50, sd = 0.1))

  fit <- npcopula(~ x + y, data = dat, neval = 4, nmulti = 1)

  expect_s3_class(fit, "npcopula")
  expect_equal(nrow(as.data.frame(fit)), 16L)
  expect_identical(fit$target, "distribution")
  expect_identical(fit$evaluation, "grid")
  expect_true(isTRUE(fit$u.auto))
  expect_equal(fit$grid.dim, c(4L, 4L))
  expect_true(inherits(fit$bws, "dbandwidth"))
  expect_true(is.data.frame(fit$eval))
  expect_equal(as.data.frame(fit), fit$eval, ignore_attr = TRUE)

  fixed.fit <- npcopula(
    ~ x + y,
    data = dat,
    bws = c(0.7, 0.7),
    bandwidth.compute = FALSE,
    neval = 3
  )
  expect_s3_class(fixed.fit, "npcopula")
  expect_equal(nrow(as.data.frame(fixed.fit)), 9L)
  expect_true(inherits(fixed.fit$bws, "dbandwidth"))

  u <- data.frame(x = c(0.25, 0.75), y = c(0.25, 0.75))
  explicit.u.fit <- npcopula(
    ~ x + y,
    data = dat,
    bws = c(0.7, 0.7),
    bandwidth.compute = FALSE,
    u = u,
    n.quasi.inv = 40
  )
  expect_s3_class(explicit.u.fit, "npcopula")
  expect_equal(nrow(as.data.frame(explicit.u.fit)), 4L)
  expect_true(inherits(explicit.u.fit$bws, "dbandwidth"))
})

test_that("npcopula formula route can request sample evaluation and density target", {
  set.seed(43)
  dat <- data.frame(x = rnorm(40), y = 0.5 * rnorm(40) + rnorm(40, sd = 0.1))

  sample.fit <- npcopula(~ x + y, data = dat, evaluation = "sample", nmulti = 1)
  expect_equal(nrow(as.data.frame(sample.fit)), nrow(dat))
  expect_identical(sample.fit$evaluation, "sample")

  dens.fit <- npcopula(~ x + y, data = dat, target = "density", neval = 3, nmulti = 1)
  expect_equal(nrow(as.data.frame(dens.fit)), 9L)
  expect_identical(dens.fit$target, "density")
  expect_true(inherits(dens.fit$bws, "bandwidth"))

  dens.fixed <- npcopula(
    ~ x + y,
    data = dat,
    target = "density",
    bws = c(0.7, 0.7),
    bandwidth.compute = FALSE,
    neval = 3
  )
  expect_equal(nrow(as.data.frame(dens.fixed)), 9L)
  expect_identical(dens.fixed$target, "density")
  expect_true(inherits(dens.fixed$bws, "bandwidth"))
})

test_that("npcopula formula route rejects unsafe automatic high-dimensional grids", {
  dat <- data.frame(x = rnorm(20), y = rnorm(20), z = rnorm(20))

  expect_error(
    npcopula(~ x + y + z, data = dat, nmulti = 1),
    "automatic copula probability grids"
  )
})

test_that("npcopula target conflicts with explicit bandwidth object fail clearly", {
  data("faithful")
  bw <- npudistbw(dat=faithful, bws=c(0.5, 5), bandwidth.compute=FALSE)
  expect_error(
    npcopula(data = faithful, bws = bw, target = "density"),
    "target.*conflicts"
  )
})

test_that("npcopula validates probability-grid and quasi-inverse controls", {
  data("faithful")
  bw <- npudistbw(dat = faithful, bws = c(0.5, 5), bandwidth.compute = FALSE)

  expect_error(
    npcopula(data = faithful, bws = bw,
             u = data.frame(eruptions = "a", waiting = "b")),
    "numeric probability"
  )
  expect_error(
    npcopula(data = faithful, bws = bw,
             u = data.frame(eruptions = 1.2, waiting = 0.5)),
    "u must lie in \\[0,1\\]"
  )
  expect_error(
    npcopula(data = faithful, bws = bw,
             u = data.frame(eruptions = 0.5)),
    "u and bws are incompatible"
  )
  expect_error(
    npcopula(data = faithful, bws = bw, u = c(0.5, 0.5),
             n.quasi.inv = 2.5),
    "integer greater than one"
  )
  expect_error(
    npcopula(data = faithful, bws = bw, u = c(0.5, 0.5),
             n.quasi.inv = NA),
    "integer greater than one"
  )
  expect_error(
    npcopula(data = faithful, bws = bw, u = c(0.5, 0.5),
             er.quasi.inv = -1),
    "non-negative finite"
  )
  expect_error(
    npcopula(data = faithful, bws = bw, u = c(0.5, 0.5),
             er.quasi.inv = NA),
    "non-negative finite"
  )
})

test_that("predict.npcopula evaluates stored bandwidths on probability grids", {
  data("faithful")
  bw <- npudistbw(dat = faithful, bws = c(0.5, 5), bandwidth.compute = FALSE)
  u <- data.frame(eruptions = seq(0.2, 0.8, length.out = 3),
                  waiting = seq(0.2, 0.8, length.out = 3))
  fit <- npcopula(data = faithful, bws = bw, u = u, n.quasi.inv = 40)

  expect_equal(predict(fit), fitted(fit))
  expect_null(attr(fit, "bws"))

  u.new <- data.frame(eruptions = c(0.25, 0.75),
                      waiting = c(0.25, 0.75))
  expected <- npcopula(bws = fit$bws, data = faithful, u = u.new,
                       n.quasi.inv = 40)
  expect_equal(predict(fit, u = u.new, n.quasi.inv = 40), fitted(expected))

  u.alias <- data.frame(u1 = c(0.25, 0.75), u2 = c(0.25, 0.75))
  expect_equal(predict(fit, newdata = u.alias, n.quasi.inv = 40),
               fitted(expected))

  pred.object <- predict(fit, u = u.new, n.quasi.inv = 40, output = "object")
  expect_s3_class(pred.object, "npcopula")
  expect_equal(as.data.frame(predict(fit, u = u.new, n.quasi.inv = 40,
                                     output = "data")),
               as.data.frame(expected),
               ignore_attr = TRUE)

  pred.se <- predict(fit, se.fit = TRUE)
  expect_equal(pred.se$fit, fitted(fit))
  expect_equal(pred.se$se.fit, se(fit))
  expect_length(pred.se$se.fit, length(fitted(fit)))
  expect_true(all(is.finite(pred.se$se.fit)))
  expect_error(predict(fit, se.fit = "yes"),
               "'se.fit' must be TRUE or FALSE", fixed = TRUE)
})

test_that("npcopula fitted and basic plot methods work", {
  data("faithful")
  bw <- npudistbw(dat=faithful, bws=c(0.5, 5), bandwidth.compute=FALSE)
  u <- data.frame(eruptions = seq(0, 1, length.out = 4),
                  waiting = seq(0, 1, length.out = 4))
  fit <- npcopula(data = faithful, bws = bw, u = u, n.quasi.inv = 40)

  expect_equal(fitted(fit), fit$copula)
  expect_silent(plot(fit, view = "fixed"))
  expect_silent(plot(fit, view = "fixed", col = NULL))
  expect_silent(plot(fit, col = "lightblue", border = "grey40",
                     theta = 20, phi = 25, view = "fixed", shade = 0.3))
  expect_silent(plot(fit, view = "image", perspective = FALSE))
  expect_silent(plot(fit, view = "contour", perspective = FALSE))
  expect_true("copula" %in%
                getFromNamespace(".np_progress_single_line_surfaces", "np")())
  expect_error(plot(fit, view = "image", renderer = "rgl"),
               "surface views")

  one.point <- npcopula(data = faithful, bws = bw, u = c(0.5, 0.5), n.quasi.inv = 40)
  expect_equal(nrow(as.data.frame(one.point)), 1L)
  expect_equal(one.point$grid.dim, c(1L, 1L))
})

test_that("npcopula contour and image plots use descriptive default titles", {
  data("faithful")
  bw <- npudistbw(dat = faithful, bws = c(0.5, 5), bandwidth.compute = FALSE)
  u <- data.frame(eruptions = seq(0, 1, length.out = 3),
                  waiting = seq(0, 1, length.out = 3))
  fit <- npcopula(data = faithful, bws = bw, u = u, n.quasi.inv = 40)

  captured <- new.env(parent = emptyenv())
  captured$contour <- NULL
  captured$image <- NULL
  pdf(file = tempfile(fileext = ".pdf"))
  on.exit(dev.off(), add = TRUE)
  trace(
    what = "contour",
    where = asNamespace("graphics"),
    tracer = bquote({
      assign("contour", as.list(match.call()), envir = .(captured))
    }),
    print = FALSE
  )
  on.exit(
    try(untrace("contour", where = asNamespace("graphics")), silent = TRUE),
    add = TRUE
  )
  trace(
    what = "image",
    where = asNamespace("graphics"),
    tracer = bquote({
      assign("image", as.list(match.call()), envir = .(captured))
    }),
    print = FALSE
  )
  on.exit(
    try(untrace("image", where = asNamespace("graphics")), silent = TRUE),
    add = TRUE
  )

  expect_silent(plot(fit, view = "contour", perspective = FALSE))
  expect_silent(plot(fit, view = "image", perspective = FALSE))
  expect_identical(eval(captured$contour$main), "Copula Contour")
  expect_identical(eval(captured$image$main), "Contour Heatmap")
})

test_that("npcopula perspective defaults follow shared surface convention", {
  plot.method <- getS3method("plot", "npcopula")
  form <- formals(plot.method)
  expect_equal(form$theta, 0.0)
  expect_equal(form$phi, 20.0)
  expect_identical(eval(form$view)[1L], "rotate")

  data("faithful")
  bw <- npudistbw(dat = faithful, bws = c(0.5, 5), bandwidth.compute = FALSE)
  u <- data.frame(eruptions = seq(0, 1, length.out = 3),
                  waiting = seq(0, 1, length.out = 3))
  fit <- npcopula(data = faithful, bws = bw, u = u, n.quasi.inv = 40)

  captured <- new.env(parent = emptyenv())
  captured$args <- list()
  pdf(file = tempfile(fileext = ".pdf"))
  on.exit(dev.off(), add = TRUE)
  trace(
    what = "persp",
    where = asNamespace("graphics"),
    tracer = bquote({
      assign("args",
             as.list(match.call()),
             envir = .(captured))
    }),
    print = FALSE
  )
  on.exit(
    try(untrace("persp", where = asNamespace("graphics")), silent = TRUE),
    add = TRUE
  )

  expect_silent(plot(fit, view = "fixed"))
  expect_equal(eval(captured$args$theta), 0.0)
  expect_equal(eval(captured$args$phi), 20.0)
  expect_identical(captured$args$ticktype, "detailed")
})

test_that("npcopula zlim caps plotted surface without changing data", {
  data("faithful")
  bw <- npudistbw(dat = faithful, bws = c(0.5, 5), bandwidth.compute = FALSE)
  u <- data.frame(eruptions = seq(0, 1, length.out = 3),
                  waiting = seq(0, 1, length.out = 3))
  fit <- npcopula(data = faithful, bws = bw, u = u, n.quasi.inv = 40)
  fit$eval$copula[1L] <- 5
  fit$copula <- fit$eval$copula

  expect_gt(max(plot(fit, output = "data")$copula), 1)
  pdf(file = tempfile(fileext = ".pdf"))
  on.exit(dev.off(), add = TRUE)
  expect_warning(plot(fit, view = "fixed", zlim = c(0, 1)), NA)
})

test_that("npcopula plot intervals work for grid surfaces", {
  set.seed(44)
  dat <- data.frame(x = rnorm(45), y = rnorm(45))
  fit <- npcopula(~ x + y, data = dat, neval = 3, nmulti = 1)

  asym <- plot(fit, output = "data", errors = "asymptotic", band = "all")
  expect_true(all(c("lower", "upper", "pointwise.lower",
                    "bonferroni.upper") %in% names(asym)))
  expect_equal(nrow(asym), nrow(as.data.frame(fit)))

  boot <- plot(fit, output = "data", errors = "bootstrap",
               bootstrap = "inid", B = 3, band = "pmzsd")
  expect_true(all(c("center", "lower", "upper") %in% names(boot)))
  expect_equal(nrow(boot), nrow(as.data.frame(fit)))

  dens.bw <- npudensbw(dat = dat, bws = c(0.7, 0.7),
                       bandwidth.compute = FALSE)
  dens.fit <- npcopula(
    bws = dens.bw,
    data = dat,
    u = data.frame(x = seq(0.2, 0.8, length.out = 3),
                   y = seq(0.2, 0.8, length.out = 3)),
    n.quasi.inv = 40
  )
  dens <- plot(dens.fit, output = "data", errors = "bootstrap",
               bootstrap = "inid", B = 3, band = "pmzsd")
  expect_true(all(c("center", "lower", "upper") %in% names(dens)))
  expect_equal(nrow(dens), nrow(as.data.frame(dens.fit)))

  sample.fit <- npcopula(~ x + y, data = dat, evaluation = "sample", nmulti = 1)
  expect_error(plot(sample.fit, errors = "asymptotic"), "grid evaluation")
})

test_that("npcopula plots mixed ordered grids on requested probability axes", {
  set.seed(42)
  n <- 80
  n.eval <- 12
  rho <- 0.99
  x <- rnorm(n)
  ylatent <- rho * x + sqrt(1 - rho^2) * rnorm(n)
  dat <- data.frame(
    x = x,
    y = ordered(as.integer(cut(
      ylatent,
      quantile(ylatent, seq(0, 1, by = 0.1)),
      include.lowest = TRUE
    )) - 1)
  )
  fit <- npcopula(~ x + y, data = dat, neval = n.eval, nmulti = 1)

  expect_equal(nrow(as.data.frame(fit)), n.eval * n.eval)
  expect_equal(length(unique(fit$u.grid$u1)), n.eval)
  expect_equal(length(unique(fit$u.grid$u2)), n.eval)
  expect_lt(length(unique(as.data.frame(fit)$u2)), n.eval)
  expect_silent(plot(fit, view = "contour", output = "data"))

  empirical <- plot(fit, view = "empirical", output = "data")
  expect_equal(nrow(empirical), n)
  expect_true(all(c("u1", "u2") %in% names(empirical)))

  all.payload <- plot(fit, view = "all", output = "data")
  expect_equal(nrow(all.payload$copula), n.eval * n.eval)
  expect_equal(nrow(all.payload$empirical), n)
  expect_equal(nrow(all.payload$density), n.eval * n.eval)
  expect_error(plot(fit, view = "all", renderer = "rgl"), "renderer='base'")

  dens <- npcopula(~ x + y, data = dat, target = "density",
                   neval = n.eval, nmulti = 1)
  expect_error(plot(dens, view = "all"),
               "requires a copula distribution object")
})
