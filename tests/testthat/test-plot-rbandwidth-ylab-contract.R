capture_rbandwidth_panel_ylabs <- function(object, xdat, ydat, ...) {
  captured <- new.env(parent = emptyenv())
  captured$ylab <- character()

  pdf(file = tempfile(fileext = ".pdf"))
  on.exit(dev.off(), add = TRUE)

  trace(
    what = "plot",
    where = asNamespace("base"),
    tracer = quote(invisible(NULL)),
    print = FALSE
  )
  on.exit(
    try(untrace("plot", where = asNamespace("base")), silent = TRUE),
    add = TRUE
  )

  trace(
    what = "plot.default",
    where = asNamespace("graphics"),
    tracer = bquote({
      assign(
        "ylab",
        c(
          get("ylab", envir = .(captured), inherits = FALSE),
          as.character(ylab)
        ),
        envir = .(captured)
      )
    }),
    print = FALSE
  )
  on.exit(
    try(untrace("plot.default", where = asNamespace("graphics")), silent = TRUE),
    add = TRUE
  )

  suppressWarnings(plot(
    object,
    xdat = xdat,
    ydat = ydat,
    perspective = FALSE,
    plot.errors.method = "none",
    plot.data.overlay = FALSE,
    plot.rug = FALSE,
    ...
  ))

  captured$ylab
}

make_regression_fixture <- function(predictors = c("g", "x"),
                                    regtype = "ll",
                                    degree = NULL,
                                    gradient.order = 1L) {
  set.seed(20260404)
  dat <- data.frame(
    y = NA_real_,
    g = factor(rep(c("a", "b"), each = 20L)),
    x = seq(0.05, 0.95, length.out = 40L)
  )
  dat$y <- 1 + 0.8 * (dat$g == "b") + sin(2 * pi * dat$x)
  formula <- stats::as.formula(paste("y ~", paste(predictors, collapse = " + ")))

  bw.args <- list(
    formula = formula,
    data = dat,
    regtype = regtype,
    bwtype = "fixed",
    bws = rep_len(c(0.4, 0.25), length(predictors)),
    bandwidth.compute = FALSE
  )
  if (!is.null(degree))
    bw.args$degree <- degree

  bw <- do.call(npregbw, bw.args)
  fit <- suppressWarnings(npreg(
    bws = bw,
    gradients = TRUE,
    gradient.order = gradient.order,
    warn.glp.gradient = FALSE
  ))

  list(
    fit = fit,
    xdat = dat[predictors],
    ydat = dat$y
  )
}

test_that("rbandwidth gradient panels use Delta for factors and d for continuous predictors", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  mixed.factor.first <- make_regression_fixture(predictors = c("g", "x"))

  mixed.labels <- capture_rbandwidth_panel_ylabs(
    mixed.factor.first$fit,
    xdat = mixed.factor.first$xdat,
    ydat = mixed.factor.first$ydat,
    gradients = TRUE,
    common.scale = FALSE
  )
  expect_identical(mixed.labels, c("Delta y / Delta g", "d y / d x"))
})

test_that("rbandwidth non-gradient default ylab is unchanged", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  fixture <- make_regression_fixture(predictors = c("x", "g"))
  labels <- capture_rbandwidth_panel_ylabs(
    fixture$fit,
    xdat = fixture$xdat,
    ydat = fixture$ydat,
    gradients = FALSE,
    common.scale = TRUE
  )

  expect_true(length(labels) >= 1L)
  expect_true(all(labels == paste("", "y")))
})

test_that("rbandwidth explicit ylab overrides remain unchanged", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  fixture <- make_regression_fixture(predictors = c("x", "g"))

  labels <- capture_rbandwidth_panel_ylabs(
    fixture$fit,
    xdat = fixture$xdat,
    ydat = fixture$ydat,
    gradients = TRUE,
    ylab = "custom",
    common.scale = TRUE
  )

  expect_true(length(labels) >= 1L)
  expect_true(all(labels == "custom"))
})

test_that("rbandwidth mixed gradient plot still returns data payloads", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)

  fixture <- make_regression_fixture(predictors = c("g", "x"))

  out <- expect_no_error(suppressWarnings(plot(
    fixture$fit,
    xdat = fixture$xdat,
    ydat = fixture$ydat,
    gradients = TRUE,
    perspective = FALSE,
    common.scale = FALSE,
    plot.behavior = "data",
    plot.errors.method = "none",
    plot.data.overlay = FALSE
  )))

  expect_type(out, "list")
  expect_length(out, 2L)
  expect_true(all(vapply(out, inherits, logical(1), "npregression")))
})
