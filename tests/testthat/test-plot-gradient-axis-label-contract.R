plot_label_contract_package <- function() {
  sub("^namespace:", "", environmentName(environment(npindexbw)))
}

skip_if_plot_label_mpi_pool_inactive <- function() {
  if (!identical(plot_label_contract_package(), "npRmpi"))
    return(invisible())
  active <- isTRUE(getOption("npRmpi.mpi.initialized", FALSE)) &&
    isTRUE(try(mpi.comm.size(1) > 1L, silent = TRUE))
  testthat::skip_if_not(active, "requires an active npRmpi session pool")
}

capture_plot_axis_labels <- function(expr) {
  pkg <- plot_label_contract_package()
  captured <- new.env(parent = emptyenv())
  captured$xlab <- list()
  captured$ylab <- list()

  grDevices::pdf(file = tempfile(fileext = ".pdf"))
  on.exit(grDevices::dev.off(), add = TRUE)
  trace(
    what = ".np_plot_merge_user_args",
    where = asNamespace(pkg),
    tracer = bquote({
      effective.args <- base.args
      if (length(user.args)) {
        user.names <- names(user.args)
        base.names <- names(effective.args)
        replaced <- intersect(
          base.names[!(is.na(base.names) | base.names == "")],
          user.names[!(is.na(user.names) | user.names == "")]
        )
        if (length(replaced))
          effective.args <- effective.args[!(base.names %in% replaced)]
        effective.args <- c(effective.args, user.args)
      }
      if (!is.null(effective.args$xlab))
        assign(
          "xlab",
          append(
            get("xlab", envir = .(captured), inherits = FALSE),
            list(effective.args$xlab)
          ),
          envir = .(captured)
        )
      if (!is.null(effective.args$ylab))
        assign(
          "ylab",
          append(
            get("ylab", envir = .(captured), inherits = FALSE),
            list(effective.args$ylab)
          ),
          envir = .(captured)
        )
    }),
    print = FALSE
  )
  on.exit(
    try(
      untrace(".np_plot_merge_user_args", where = asNamespace(pkg)),
      silent = TRUE
    ),
    add = TRUE
  )

  eval.parent(substitute(expr))
  list(xlab = unlist(captured$xlab), ylab = unlist(captured$ylab))
}

test_that("shared plot derivative labels distinguish derivatives and contrasts", {
  pkg <- plot_label_contract_package()
  label <- getFromNamespace(".np_plot_gradient_axis_label", pkg)
  conditional.label <- getFromNamespace(
    ".np_plot_conditional_gradient_axis_label",
    pkg
  )

  expect_identical(label("y", "x"), "d y / d x")
  expect_identical(label("y", "x", order = 3L), "d^3 y / d x^3")
  expect_identical(
    label("y", "g", categorical = TRUE),
    "Delta y / Delta g"
  )
  expect_identical(
    conditional.label(
      "Conditional Density",
      "x2",
      component = 3L,
      continuous = c(TRUE, FALSE, TRUE),
      gradient.order = c(2L, 3L)
    ),
    "d^3 Conditional Density / d x2^3"
  )
  expect_identical(
    conditional.label(
      "Conditional Density",
      "g",
      component = 2L,
      continuous = c(TRUE, FALSE, TRUE),
      gradient.order = c(2L, 3L)
    ),
    "Delta Conditional Density / Delta g"
  )
  expect_error(label("y", "x", order = 1.5), "positive integer")
})

test_that("npindex plot labels use the index and named marginal effects", {
  skip_if_plot_label_mpi_pool_inactive()
  old.options <- options(np.messages = FALSE)
  on.exit(options(old.options), add = TRUE)

  set.seed(20260713)
  dat <- data.frame(x1 = runif(32L, -1, 1), x2 = runif(32L, -1, 1))
  dat$y <- dat$x1 - dat$x2 + rnorm(32L, sd = 0.05)
  bw <- npindexbw(
    y ~ x1 + x2,
    data = dat,
    bws = c(1, -1, 0.4),
    bandwidth.compute = FALSE,
    regtype = "lp",
    degree = 2L
  )
  fit <- npindex(bws = bw, gradients = TRUE)

  labels <- capture_plot_axis_labels(plot(
    fit,
    gradients = TRUE,
    errors = "none",
    data_rug = FALSE
  ))
  expect_true("Index, X' beta" %in% labels$xlab)
  expect_true(all(c("d y / d x1", "d y / d x2") %in% labels$ylab))

  override <- capture_plot_axis_labels(plot(
    fit,
    gradients = TRUE,
    errors = "none",
    xlab = "custom x",
    ylab = "custom y",
    data_rug = FALSE
  ))
  expect_true(length(override$xlab) > 0L)
  expect_true(length(override$ylab) > 0L)
  expect_true(all(override$xlab == "custom x"))
  expect_true(all(override$ylab == "custom y"))
})

test_that("semiparametric plot labels are mathematical and reachable", {
  skip_if_plot_label_mpi_pool_inactive()
  old.options <- options(np.messages = FALSE)
  on.exit(options(old.options), add = TRUE)

  set.seed(20260714)
  dat <- data.frame(x = runif(32L), z = runif(32L))
  dat$y <- (1 + dat$x^2) * dat$z + rnorm(32L, sd = 0.05)
  sc.bw <- npscoefbw(
    y ~ x | z,
    data = dat,
    bws = 0.4,
    bandwidth.compute = FALSE,
    regtype = "lp",
    degree = 1L
  )
  sc.fit <- npscoef(bws = sc.bw, gradients = TRUE, errors = FALSE)
  sc.labels <- capture_plot_axis_labels(plot(
    sc.fit,
    gradients = TRUE,
    coef_index = 1L,
    perspective = FALSE,
    common_scale = FALSE,
    errors = "none",
    data_overlay = FALSE,
    data_rug = FALSE,
    neval = 7L
  ))
  expect_true("d y / d x" %in% sc.labels$ylab)

  pl.bw <- npplregbw(
    y ~ z | x,
    data = dat,
    bws = matrix(c(0.4, 0.4), nrow = 2L),
    bandwidth.compute = FALSE,
    regtype = "lp",
    degree = 1L
  )
  pl.fit <- npplreg(bws = pl.bw)
  pl.labels <- capture_plot_axis_labels(plot(
    pl.fit,
    perspective = FALSE,
    common_scale = FALSE,
    errors = "none",
    data_overlay = FALSE,
    data_rug = FALSE,
    neval = 7L
  ))
  expect_true(length(pl.labels$ylab) > 0L)
  expect_true(all(pl.labels$ylab == "y"))
  expect_error(
    plot(pl.fit, gradients = TRUE),
    "gradients not supported with partially linear models"
  )
})

test_that("conditional plot labels map components by type and order", {
  skip_if_plot_label_mpi_pool_inactive()
  old.options <- options(np.messages = FALSE)
  on.exit(options(old.options), add = TRUE)

  set.seed(20260715)
  dat <- data.frame(
    x = seq(0.02, 0.98, length.out = 32L),
    g = factor(rep(c("a", "b"), length.out = 32L))
  )
  dat$y <- sin(2 * pi * dat$x) + 0.4 * (dat$g == "b") +
    rnorm(32L, sd = 0.08)
  bw <- npcdensbw(
    y ~ x + g,
    data = dat,
    bws = c(0.4, 0.25, 0.4),
    bandwidth.compute = FALSE,
    bwtype = "fixed"
  )
  labels <- capture_plot_axis_labels(plot(
    bw,
    gradients = TRUE,
    perspective = FALSE,
    common_scale = FALSE,
    errors = "none",
    neval = 7L
  ))
  expect_true(all(c(
    "d Conditional Density / d x",
    "Delta Conditional Density / Delta g"
  ) %in% labels$ylab))

  lp.dat <- data.frame(x = seq(-1, 1, length.out = 32L))
  lp.dat$y <- lp.dat$x^3 + rnorm(32L, sd = 0.03)
  lp.bw <- npcdensbw(
    y ~ x,
    data = lp.dat,
    bws = c(0.5, 0.5),
    bandwidth.compute = FALSE,
    bwtype = "fixed",
    regtype = "lp",
    degree = 2L
  )
  lp.labels <- capture_plot_axis_labels(plot(
    lp.bw,
    gradients = TRUE,
    gradient_order = 2L,
    perspective = FALSE,
    errors = "none",
    neval = 7L
  ))
  expect_true("d^2 Conditional Density / d x^2" %in% lp.labels$ylab)
})
