with_nprmpi_bindings <- function(bindings, code) {
  code <- substitute(code)
  ns <- asNamespace("npRmpi")
  old <- lapply(names(bindings), function(name) get(name, envir = ns, inherits = FALSE))
  names(old) <- names(bindings)

  for (name in names(bindings)) {
    was_locked <- bindingIsLocked(name, ns)
    if (was_locked) {
      unlockBinding(name, ns)
    }
    assign(name, bindings[[name]], envir = ns)
    if (was_locked) {
      lockBinding(name, ns)
    }
  }

  on.exit({
    for (name in names(old)) {
      was_locked <- bindingIsLocked(name, ns)
      if (was_locked) {
        unlockBinding(name, ns)
      }
      assign(name, old[[name]], envir = ns)
      if (was_locked) {
        lockBinding(name, ns)
      }
    }
  }, add = TRUE)

  eval(code, envir = parent.frame())
}

test_that("regression bootstrap target labels format compactly", {
  fmt <- getFromNamespace(".np_plot_regression_bootstrap_target_label", "npRmpi")

  bws <- list(
    ndim = 2L,
    xnames = c("x1", "x2")
  )

  expect_identical(fmt(bws = bws, slice.index = 0L, gradients = FALSE), "surf 1/2")
  expect_identical(fmt(bws = bws, slice.index = 1L, gradients = FALSE), "x1 1/2")
  expect_identical(fmt(bws = bws, slice.index = 2L, gradients = TRUE), "grad x2 2/2")
})

test_that("conditional bootstrap target labels format compactly", {
  fmt <- getFromNamespace(".np_plot_conditional_bootstrap_target_label", "npRmpi")

  bws <- list(
    xndim = 2L,
    yndim = 1L,
    xnames = c("x1", "x2"),
    ynames = "y"
  )

  expect_identical(fmt(bws = bws, slice.index = 0L, gradients = FALSE), "surf 1/3")
  expect_identical(fmt(bws = bws, slice.index = 1L, gradients = FALSE), "x1 1/3")
  expect_identical(fmt(bws = bws, slice.index = 3L, gradients = FALSE), "y 3/3")
  expect_identical(fmt(bws = bws, slice.index = 2L, gradients = TRUE, gradient.index = 1L), "grad x1 on x2 2/3")
})

test_that("single-index bootstrap target labels format compactly", {
  fmt <- getFromNamespace(".np_plot_singleindex_bootstrap_target_label", "npRmpi")

  expect_null(fmt(gradients = FALSE))
  expect_null(fmt(gradients = TRUE))
})

test_that("smooth coefficient bootstrap target labels format compactly", {
  fmt <- getFromNamespace(".np_plot_scoef_bootstrap_target_label", "npRmpi")

  bws <- list(
    xndim = 1L,
    zndim = 1L,
    xnames = "x",
    znames = "z"
  )

  expect_identical(fmt(bws = bws, slice.index = 0L), "surf 1/2")
  expect_identical(fmt(bws = bws, slice.index = 1L), "x 1/2")
  expect_identical(fmt(bws = bws, slice.index = 2L), "z 2/2")
})

test_that("regression helper labels carry target context for block bootstrap phases", {
  fn <- getFromNamespace("compute.bootstrap.errors.rbandwidth", "npRmpi")
  captured <- new.env(parent = emptyenv())
  captured$prep <- character()
  captured$progress <- character()
  captured$interval <- character()

  res <- with_nprmpi_bindings(
    list(
      .np_plot_activity_begin = function(label) {
        captured$prep <- c(captured$prep, label)
        list(label = label)
      },
      .np_plot_activity_end = function(state) invisible(NULL),
      .np_inid_boot_from_regression = function(..., prep.label = NULL, progress.label = NULL) {
        captured$progress <- c(captured$progress, prep.label, progress.label)
        list(t = matrix(0, nrow = 2L, ncol = 3L), t0 = c(0, 0, 0))
      },
      .np_plot_bootstrap_interval_summary = function(boot.t, t0, alpha, band.type, progress.label = NULL) {
        captured$interval <- c(captured$interval, progress.label)
        list(err = matrix(0, nrow = ncol(boot.t), ncol = 2L), all.err = list())
      },
      .npRmpi_profile_bootstrap_begin = function(...) list(),
      .npRmpi_profile_finalize_bootstrap = function(boot.err, bxp, boot.all.err, ctx) {
        list(boot.err = boot.err, bxp = bxp, boot.all.err = boot.all.err)
      }
    ),
    fn(
      xdat = data.frame(x1 = c(0, 1, 2), x2 = c(1, 2, 3)),
      ydat = c(1, 2, 3),
      exdat = data.frame(x1 = c(0, 1, 2), x2 = c(1, 2, 3)),
      gradients = TRUE,
      gradient.order = 1L,
      slice.index = 2L,
      plot.errors.boot.method = "geom",
      plot.errors.boot.blocklen = 2L,
      plot.errors.boot.num = 2L,
      plot.errors.center = "estimate",
      plot.errors.type = "all",
      plot.errors.alpha = 0.05,
      progress.target = "grad x2 2/2",
      bws = list(
        type = "fixed",
        ndim = 2L,
        xnames = c("x1", "x2"),
        xdati = list(icon = c(TRUE, TRUE), iord = c(FALSE, FALSE), iuno = c(FALSE, FALSE))
      )
    )
  )

  expect_true(is.list(res))
  expect_true(any(grepl("Preparing plot bootstrap geom \\(grad x2 2/2\\)", captured$prep)))
  expect_true(any(grepl("Preparing plot bootstrap geom \\(grad x2 2/2\\)", captured$progress)))
  expect_true(any(grepl("Plot bootstrap \\(grad x2 2/2\\)", captured$progress)))
  expect_true(any(grepl("Constructing bootstrap all bands \\(grad x2 2/2\\)", captured$interval)))
})

test_that("conditional helper labels carry target context for block bootstrap phases", {
  fn <- getFromNamespace("compute.bootstrap.errors.conbandwidth", "npRmpi")
  captured <- new.env(parent = emptyenv())
  captured$prep <- character()
  captured$progress <- character()
  captured$interval <- character()

  res <- with_nprmpi_bindings(
    list(
      .np_plot_activity_begin = function(label) {
        captured$prep <- c(captured$prep, label)
        list(label = label)
      },
      .np_plot_activity_end = function(state) invisible(NULL),
      .np_inid_boot_from_conditional_gradient_local = function(..., progress.label = NULL) {
        captured$progress <- c(captured$progress, progress.label)
        list(t = matrix(0, nrow = 2L, ncol = 3L), t0 = c(0, 0, 0))
      },
      .np_plot_bootstrap_interval_summary = function(boot.t, t0, alpha, band.type, progress.label = NULL) {
        captured$interval <- c(captured$interval, progress.label)
        list(err = matrix(0, nrow = ncol(boot.t), ncol = 2L), all.err = list())
      },
      .npRmpi_profile_bootstrap_begin = function(...) list(),
      .npRmpi_profile_finalize_bootstrap = function(boot.err, bxp, boot.all.err, ctx) {
        list(boot.err = boot.err, bxp = bxp, boot.all.err = boot.all.err)
      },
      .npRmpi_with_local_bootstrap = function(expr) expr
    ),
    fn(
      xdat = data.frame(x1 = c(0, 1, 2), x2 = c(1, 2, 3)),
      ydat = data.frame(y = c(0.1, 0.2, 0.3)),
      exdat = data.frame(x1 = c(0, 1, 2), x2 = c(1, 2, 3)),
      eydat = data.frame(y = c(0.1, 0.2, 0.3)),
      cdf = FALSE,
      quantreg = FALSE,
      tau = NULL,
      gradients = TRUE,
      gradient.index = 1L,
      slice.index = 2L,
      plot.errors.boot.method = "geom",
      plot.errors.boot.nonfixed = "exact",
      plot.errors.boot.blocklen = 2L,
      plot.errors.boot.num = 2L,
      plot.errors.center = "estimate",
      plot.errors.type = "all",
      plot.errors.alpha = 0.05,
      progress.target = "grad x1 on x2 2/3",
      bws = list(
        type = "fixed",
        xndim = 2L,
        yndim = 1L,
        xnames = c("x1", "x2"),
        ynames = "y",
        xdati = list(),
        ydati = list()
      )
    )
  )

  expect_true(is.list(res))
  expect_true(any(grepl("Preparing plot bootstrap geom \\(grad x1 on x2 2/3\\)", captured$prep)))
  expect_true(any(grepl("Plot bootstrap \\(grad x1 on x2 2/3\\)", captured$progress)))
  expect_true(any(grepl("Constructing bootstrap all bands \\(grad x1 on x2 2/3\\)", captured$interval)))
})

test_that("conditional exact ksum wrapper preserves user-facing progress label", {
  fn <- getFromNamespace(".np_inid_boot_from_ksum_conditional", "npRmpi")
  captured <- new.env(parent = emptyenv())
  captured$progress <- character()

  res <- with_nprmpi_bindings(
    list(
      .np_inid_boot_from_ksum_conditional_exact = function(..., progress.label = NULL) {
        captured$progress <- c(captured$progress, progress.label)
        list(t = matrix(0, nrow = 2L, ncol = 2L), t0 = c(0, 0))
      }
    ),
    fn(
      xdat = data.frame(x = c(0, 1, 2)),
      ydat = data.frame(y = c(0.1, 0.2, 0.3)),
      exdat = data.frame(x = c(0, 1)),
      eydat = data.frame(y = c(0.1, 0.2)),
      bws = list(type = "adaptive_nn"),
      B = 2L,
      cdf = FALSE,
      counts.drawer = function(start, stopi) matrix(1L, nrow = 3L, ncol = stopi - start + 1L),
      progress.label = "Plot bootstrap"
    )
  )

  expect_true(is.list(res))
  expect_identical(captured$progress, "Plot bootstrap")
})

test_that("conditional bootstrap exact path forwards target label from compute helper", {
  fn <- getFromNamespace("compute.bootstrap.errors.conbandwidth", "npRmpi")
  captured <- new.env(parent = emptyenv())
  captured$prep <- character()
  captured$progress <- character()
  captured$interval <- character()

  res <- with_nprmpi_bindings(
    list(
      .np_plot_activity_begin = function(label) {
        captured$prep <- c(captured$prep, label)
        list(label = label)
      },
      .np_plot_activity_end = function(state) invisible(NULL),
      .np_inid_boot_from_ksum_conditional = function(..., progress.label = NULL) {
        captured$progress <- c(captured$progress, progress.label)
        list(t = matrix(0, nrow = 2L, ncol = 3L), t0 = c(0, 0, 0))
      },
      .np_plot_bootstrap_interval_summary = function(boot.t, t0, alpha, band.type, progress.label = NULL) {
        captured$interval <- c(captured$interval, progress.label)
        list(err = matrix(0, nrow = ncol(boot.t), ncol = 2L), all.err = list())
      },
      .npRmpi_profile_bootstrap_begin = function(...) list(),
      .npRmpi_profile_finalize_bootstrap = function(boot.err, bxp, boot.all.err, ctx) {
        list(boot.err = boot.err, bxp = bxp, boot.all.err = boot.all.err)
      },
      .npRmpi_with_local_bootstrap = function(expr) expr,
      .np_con_inid_ksum_eligible = function(bws) TRUE
    ),
    fn(
      xdat = data.frame(x = c(0, 1, 2)),
      ydat = data.frame(y = c(0.1, 0.2, 0.3)),
      exdat = data.frame(x = c(0, 1, 2)),
      eydat = data.frame(y = c(0.1, 0.2, 0.3)),
      cdf = FALSE,
      quantreg = FALSE,
      tau = NULL,
      gradients = FALSE,
      gradient.index = NULL,
      slice.index = 0L,
      plot.errors.boot.method = "geom",
      plot.errors.boot.nonfixed = "exact",
      plot.errors.boot.blocklen = 2L,
      plot.errors.boot.num = 2L,
      plot.errors.center = "estimate",
      plot.errors.type = "all",
      plot.errors.alpha = 0.05,
      progress.target = NULL,
      bws = list(
        type = "adaptive_nn",
        xndim = 1L,
        yndim = 1L,
        xnames = "x",
        ynames = "y",
        xdati = list(),
        ydati = list()
      )
    )
  )

  expect_true(is.list(res))
  expect_true(any(grepl("^Preparing plot bootstrap geom$", captured$prep)))
  expect_identical(captured$progress, "Plot bootstrap")
  expect_true(any(grepl("^Constructing bootstrap all bands$", captured$interval)))
})

test_that("single-index helper labels carry target context for bootstrap phases", {
  fn <- getFromNamespace("compute.bootstrap.errors.sibandwidth", "npRmpi")
  captured <- new.env(parent = emptyenv())
  captured$prep <- character()
  captured$progress <- character()
  captured$interval <- character()

  res <- with_nprmpi_bindings(
    list(
      .np_plot_activity_begin = function(label) {
        captured$prep <- c(captured$prep, label)
        list(label = label)
      },
      .np_plot_activity_end = function(state) invisible(NULL),
      .np_inid_boot_from_index = function(..., progress.label = NULL) {
        captured$progress <- c(captured$progress, progress.label)
        list(t = matrix(0, nrow = 2L, ncol = 3L), t0 = c(0, 0, 0))
      },
      .np_plot_bootstrap_interval_summary = function(boot.t, t0, alpha, band.type, progress.label = NULL) {
        captured$interval <- c(captured$interval, progress.label)
        list(err = matrix(0, nrow = ncol(boot.t), ncol = 2L), all.err = list())
      },
      .npRmpi_profile_bootstrap_begin = function(...) list(),
      .npRmpi_profile_finalize_bootstrap = function(boot.err, bxp, boot.all.err, ctx) {
        list(boot.err = boot.err, bxp = bxp, boot.all.err = boot.all.err)
      },
      .npRmpi_with_local_bootstrap = function(expr) expr
    ),
    fn(
      xdat = data.frame(x1 = c(0, 1, 2), x2 = c(1, 2, 3)),
      ydat = c(0.1, 0.2, 0.3),
      gradients = FALSE,
      plot.errors.boot.method = "geom",
      plot.errors.boot.blocklen = 2L,
      plot.errors.boot.num = 2L,
      plot.errors.center = "estimate",
      plot.errors.type = "all",
      plot.errors.alpha = 0.05,
      bws = list(type = "fixed", beta = c(1, 1))
    )
  )

  expect_true(is.list(res))
  expect_true(any(grepl("^Preparing plot bootstrap geom$", captured$prep)))
  expect_true(any(grepl("^Plot bootstrap$", captured$progress)))
  expect_true(any(grepl("^Constructing bootstrap all bands$", captured$interval)))
})

test_that("unconditional helper keeps user-facing bootstrap labels in MPI fanout", {
  fn <- getFromNamespace(".np_inid_boot_from_ksum_unconditional", "npRmpi")
  captured <- new.env(parent = emptyenv())
  captured$labels <- character()

  res <- with_nprmpi_bindings(
    list(
      .np_plot_kernel_weights_direct = function(...) matrix(
        c(1, 2, 3, 4, 5, 6),
        nrow = 3L,
        ncol = 2L
      ),
      .npRmpi_bootstrap_tune_chunk_size = function(...) 2L,
      .npRmpi_bootstrap_fanout_enabled = function(...) TRUE,
      .npRmpi_bootstrap_chunk_tasks = function(B, chunk.size) list(list(start = 1L, bsz = as.integer(B), seed = 1L)),
      .npRmpi_bootstrap_run_fanout = function(tasks,
                                              worker,
                                              ncol.out,
                                              what = "bootstrap",
                                              progress.label = NULL,
                                              ...) {
        captured$labels <- c(captured$labels, progress.label)
        matrix(1, nrow = as.integer(tasks[[1]]$bsz), ncol = ncol.out)
      },
      .npRmpi_bootstrap_fail_or_fallback = function(...) stop("unexpected fallback")
    ),
    fn(
      xdat = data.frame(y = c(0.1, 0.2, 0.3)),
      exdat = data.frame(y = c(0.1, 0.2)),
      bws = list(type = "fixed"),
      B = 2L,
      operator = "normal",
      counts.drawer = function(start, stop) seq.int(start, stop)
    )
  )

  expect_true(is.list(res))
  expect_true(any(captured$labels == "Plot bootstrap block"))
  expect_false(any(grepl("inid-ksum-unconditional", captured$labels, fixed = TRUE)))
})

test_that("smooth coefficient helper labels carry target context for bootstrap phases", {
  fn <- getFromNamespace("compute.bootstrap.errors.scbandwidth", "npRmpi")
  captured <- new.env(parent = emptyenv())
  captured$prep <- character()
  captured$progress <- character()
  captured$interval <- character()

  res <- with_nprmpi_bindings(
    list(
      .np_plot_activity_begin = function(label) {
        captured$prep <- c(captured$prep, label)
        list(label = label)
      },
      .np_plot_activity_end = function(state) invisible(NULL),
      .np_inid_boot_from_scoef = function(..., progress.label = NULL) {
        captured$progress <- c(captured$progress, progress.label)
        list(t = matrix(0, nrow = 2L, ncol = 3L), t0 = c(0, 0, 0))
      },
      .np_plot_bootstrap_interval_summary = function(boot.t, t0, alpha, band.type, progress.label = NULL) {
        captured$interval <- c(captured$interval, progress.label)
        list(err = matrix(0, nrow = ncol(boot.t), ncol = 2L), all.err = list())
      },
      .npRmpi_profile_bootstrap_begin = function(...) list(),
      .npRmpi_profile_finalize_bootstrap = function(boot.err, bxp, boot.all.err, ctx) {
        list(boot.err = boot.err, bxp = bxp, boot.all.err = boot.all.err)
      },
      .npRmpi_with_local_bootstrap = function(expr) expr
    ),
    fn(
      xdat = data.frame(x = c(0, 1, 2)),
      ydat = c(0.1, 0.2, 0.3),
      zdat = data.frame(z = c(1, 2, 3)),
      exdat = data.frame(x = c(0, 1, 2)),
      ezdat = data.frame(z = c(1, 2, 3)),
      gradients = FALSE,
      slice.index = 2L,
      progress.target = "z 2/2",
      plot.errors.boot.method = "geom",
      plot.errors.boot.blocklen = 2L,
      plot.errors.boot.num = 2L,
      plot.errors.center = "estimate",
      plot.errors.type = "all",
      plot.errors.alpha = 0.05,
      bws = list(
        xndim = 1L,
        zndim = 1L,
        xnames = "x",
        znames = "z",
        xdati = list(iord = FALSE, iuno = FALSE),
        zdati = list(iord = FALSE, iuno = FALSE)
      )
    )
  )

  expect_true(is.list(res))
  expect_true(any(grepl("Preparing plot bootstrap geom \\(z 2/2\\)", captured$prep)))
  expect_true(any(grepl("Plot bootstrap \\(z 2/2\\)", captured$progress)))
  expect_true(any(grepl("Constructing bootstrap all bands \\(z 2/2\\)", captured$interval)))
})

test_that("partially linear helper labels carry target context for bootstrap phases", {
  fn <- getFromNamespace("compute.bootstrap.errors.plbandwidth", "npRmpi")
  captured <- new.env(parent = emptyenv())
  captured$prep <- character()
  captured$progress <- character()
  captured$interval <- character()

  res <- with_nprmpi_bindings(
    list(
      .np_plot_activity_begin = function(label) {
        captured$prep <- c(captured$prep, label)
        list(label = label)
      },
      .np_plot_activity_end = function(state) invisible(NULL),
      .np_inid_boot_from_plreg = function(..., progress.label = NULL) {
        captured$progress <- c(captured$progress, progress.label)
        list(t = matrix(0, nrow = 2L, ncol = 3L), t0 = c(0, 0, 0))
      },
      .np_plot_bootstrap_interval_summary = function(boot.t, t0, alpha, band.type, progress.label = NULL) {
        captured$interval <- c(captured$interval, progress.label)
        list(err = matrix(0, nrow = ncol(boot.t), ncol = 2L), all.err = list())
      },
      .npRmpi_profile_bootstrap_begin = function(...) list(),
      .npRmpi_profile_finalize_bootstrap = function(boot.err, bxp, boot.all.err, ctx) {
        list(boot.err = boot.err, bxp = bxp, boot.all.err = boot.all.err)
      }
    ),
    fn(
      xdat = data.frame(x = c(0, 1, 2)),
      ydat = c(0.1, 0.2, 0.3),
      zdat = data.frame(z = c(1, 2, 3)),
      exdat = data.frame(x = c(0, 1, 2)),
      ezdat = data.frame(z = c(1, 2, 3)),
      gradients = FALSE,
      slice.index = 2L,
      progress.target = "z 2/2",
      plot.errors.boot.method = "geom",
      plot.errors.boot.blocklen = 2L,
      plot.errors.boot.num = 2L,
      plot.errors.center = "estimate",
      plot.errors.type = "all",
      plot.errors.alpha = 0.05,
      bws = list(
        xndim = 1L,
        zndim = 1L,
        xnames = "x",
        znames = "z",
        xdati = list(iord = FALSE, iuno = FALSE),
        zdati = list(iord = FALSE, iuno = FALSE)
      )
    )
  )

  expect_true(is.list(res))
  expect_true(any(grepl("Preparing plot bootstrap geom \\(z 2/2\\)", captured$prep)))
  expect_true(any(grepl("Plot bootstrap \\(z 2/2\\)", captured$progress)))
  expect_true(any(grepl("Constructing bootstrap all bands \\(z 2/2\\)", captured$interval)))
})
