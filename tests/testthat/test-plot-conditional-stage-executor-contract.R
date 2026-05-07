test_that("conditional plot bootstrap executor preserves staged density/distribution contracts", {
  compute <- getFromNamespace("compute.bootstrap.errors.conbandwidth", "npRmpi")
  called <- new.env(parent = emptyenv())
  called$bootstrap <- list()
  called$proper <- list()

  make_bws <- function(regtype.engine = "lp") {
    structure(
      list(
        type = "fixed",
        xndim = 1L,
        yndim = 1L,
        regtype.engine = regtype.engine,
        cxkertype = "gaussian",
        cykertype = "gaussian",
        cxkerorder = 2L,
        cykerorder = 2L,
        uxkertype = "aitchisonaitken",
        uykertype = "aitchisonaitken",
        oxkertype = "wangvanryzin",
        oykertype = "wangvanryzin",
        xdati = list(all.lev = list("x"), all.ulev = list(1), iord = FALSE),
        ydati = list(all.lev = list("y"), all.ulev = list(1), iord = FALSE)
      ),
      class = "conbandwidth"
    )
  }

  local_mocked_bindings(
    .np_plot_activity_begin = function(label) list(label = label),
    .np_plot_activity_end = function(activity) invisible(TRUE),
    .np_plot_stage_progress_begin = function(total, label) list(total = total, label = label),
    .np_plot_progress_tick = function(state, done, force = FALSE) state,
    .np_plot_progress_end = function(state) invisible(TRUE),
    .np_plot_require_bws = function(bws, where) invisible(TRUE),
    .npRmpi_profile_bootstrap_begin = function(...) NULL,
    .npRmpi_profile_finalize_bootstrap = function(boot.err, bxp, boot.all.err, ctx) {
      list(boot.err = boot.err, bxp = bxp, boot.all.err = boot.all.err)
    },
    .npRmpi_with_local_bootstrap = function(expr) force(expr),
    .np_plot_boot_factor_boxplots = function(...) NULL,
    .np_con_inid_ksum_eligible = function(bws) TRUE,
    .np_plot_inid_fastpath_enabled = function() TRUE,
    .npRmpi_plot_inid_ksum_fastpath_enabled = function() TRUE,
    .np_inid_boot_from_ksum_conditional = function(xdat,
                                                   ydat,
                                                   exdat,
                                                   eydat,
                                                   bws,
                                                   B,
                                                   cdf,
                                                   counts = NULL,
                                                   counts.drawer = NULL,
                                                   progress.label = NULL) {
      called$bootstrap[[length(called$bootstrap) + 1L]] <- list(
        family = if (isTRUE(cdf)) "npcdist" else "npcdens",
        regtype.engine = bws$regtype.engine,
        B = as.integer(B),
        neval = nrow(as.data.frame(exdat)),
        progress.label = progress.label
      )
      t0 <- seq_len(nrow(as.data.frame(exdat)))
      t <- matrix(
        rep(t0, each = as.integer(B)) + seq_len(as.integer(B)) / 100,
        nrow = as.integer(B),
        ncol = length(t0)
      )
      list(t = t, t0 = t0)
    },
    .np_condens_prepare_proper_plan = function(object, proper.control) {
      list(supported = TRUE, reason = "ok", family = "npcdens")
    },
    .np_condist_prepare_proper_plan = function(object, proper.control) {
      list(supported = TRUE, reason = "ok", family = "npcdist")
    },
    .npRmpi_project_conditional_proper_values = function(values,
                                                         plan,
                                                         cdf,
                                                         progress.label = NULL,
                                                         what = NULL) {
      called$proper[[length(called$proper) + 1L]] <- list(
        family = if (isTRUE(cdf)) "npcdist" else "npcdens",
        vector = is.null(dim(values)),
        rows = if (is.null(dim(values))) length(values) else nrow(values),
        cols = if (is.null(dim(values))) 1L else ncol(values),
        progress.label = progress.label
      )
      values
    },
    .package = "npRmpi"
  )

  xdat <- data.frame(x = seq(0.1, 0.9, length.out = 6L))
  ydat <- data.frame(y = seq(0.2, 1.2, length.out = 6L))
  exdat <- data.frame(x = c(0.25, 0.50, 0.75))
  eydat <- data.frame(y = c(0.30, 0.60, 0.90))

  dens <- compute(
    xdat = xdat,
    ydat = ydat,
    exdat = exdat,
    eydat = eydat,
    cdf = FALSE,
    quantreg = FALSE,
    tau = 0.5,
    gradients = FALSE,
    gradient.index = 0L,
    slice.index = 1L,
    plot.errors.boot.method = "inid",
    plot.errors.boot.nonfixed = "exact",
    plot.errors.boot.blocklen = NULL,
    plot.errors.boot.num = 5L,
    plot.errors.center = "estimate",
    plot.errors.type = "pmzsd",
    plot.errors.alpha = 0.05,
    proper = TRUE,
    proper.method = "project",
    proper.control = list(),
    bws = make_bws("lp")
  )
  dist <- compute(
    xdat = xdat,
    ydat = ydat,
    exdat = exdat,
    eydat = eydat,
    cdf = TRUE,
    quantreg = FALSE,
    tau = 0.5,
    gradients = FALSE,
    gradient.index = 0L,
    slice.index = 1L,
    plot.errors.boot.method = "inid",
    plot.errors.boot.nonfixed = "exact",
    plot.errors.boot.blocklen = NULL,
    plot.errors.boot.num = 5L,
    plot.errors.center = "estimate",
    plot.errors.type = "pmzsd",
    plot.errors.alpha = 0.05,
    proper = TRUE,
    proper.method = "isotonic",
    proper.control = list(),
    bws = make_bws("lp")
  )

  expect_named(dens, c("boot.err", "bxp", "boot.all.err"))
  expect_named(dist, c("boot.err", "bxp", "boot.all.err"))
  expect_equal(dim(dens$boot.err), c(3L, 3L))
  expect_equal(dim(dist$boot.err), c(3L, 3L))

  expect_length(called$bootstrap, 2L)
  expect_identical(vapply(called$bootstrap, `[[`, character(1), "family"),
                   c("npcdens", "npcdist"))
  expect_true(all(vapply(called$bootstrap, `[[`, character(1), "regtype.engine") == "lp"))
  expect_true(all(vapply(called$bootstrap, `[[`, integer(1), "B") == 5L))

  expect_length(called$proper, 4L)
  expect_identical(vapply(called$proper, `[[`, character(1), "family"),
                   c("npcdens", "npcdens", "npcdist", "npcdist"))
  expect_identical(vapply(called$proper, `[[`, logical(1), "vector"),
                   c(TRUE, FALSE, TRUE, FALSE))
})
