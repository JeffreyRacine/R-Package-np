capture_rbandwidth_panel_ylabs <- function(object, xdat, ydat, ...) {
  captured <- new.env(parent = emptyenv())
  captured$ylab <- character()

  pdf(file = tempfile(fileext = ".pdf"))
  on.exit(dev.off(), add = TRUE)

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
    errors = "none",
    data_overlay = FALSE,
    data_rug = FALSE,
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
    txdat = dat[predictors],
    tydat = dat$y,
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

test_that("rbandwidth factor gradient panels use Delta labels", {
  mixed.factor.first <- make_regression_fixture(predictors = c("g", "x"))

  mixed.labels <- capture_rbandwidth_panel_ylabs(
    mixed.factor.first$fit,
    xdat = mixed.factor.first$xdat,
    ydat = mixed.factor.first$ydat,
    gradients = TRUE,
    common_scale = FALSE
  )
  expect_true("Delta y / Delta g" %in% mixed.labels)
})

test_that("rbandwidth non-gradient default ylab is unchanged", {
  fixture <- make_regression_fixture(predictors = c("x", "g"))
  labels <- capture_rbandwidth_panel_ylabs(
    fixture$fit,
    xdat = fixture$xdat,
    ydat = fixture$ydat,
    gradients = FALSE,
    common_scale = TRUE
  )

  expect_true(length(labels) >= 1L)
  expect_true(all(labels == paste("", "y")))
})

test_that("rbandwidth explicit ylab overrides remain unchanged", {
  fixture <- make_regression_fixture(predictors = c("x", "g"))

  labels <- capture_rbandwidth_panel_ylabs(
    fixture$fit,
    xdat = fixture$xdat,
    ydat = fixture$ydat,
    gradients = TRUE,
    ylab = "custom",
    common_scale = TRUE
  )

  expect_true(length(labels) >= 1L)
  expect_true(all(labels == "custom"))
})

test_that("rbandwidth mixed gradient plot still returns data payloads", {
  fixture <- make_regression_fixture(predictors = c("g", "x"))

  out <- expect_no_error(suppressWarnings(plot(
    fixture$fit,
    xdat = fixture$xdat,
    ydat = fixture$ydat,
    gradients = TRUE,
    perspective = FALSE,
    common_scale = FALSE,
    output = "data",
    errors = "none",
    data_overlay = FALSE
  )))

  expect_type(out, "list")
  expect_length(out, 2L)
  expect_true(all(vapply(out, inherits, logical(1), "npregression")))
})
