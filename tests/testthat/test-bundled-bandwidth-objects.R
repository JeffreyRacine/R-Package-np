test_that("bundled bandwidth objects retain reusable canonical metadata", {
  pkg <- "npRmpi"
  skip_on_cran()
  if (!spawn_mpi_slaves(1L)) skip("MPI slaves unavailable")
  on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  withr::local_options(list(np.messages = FALSE, np.tree = FALSE))
  withr::local_seed(123)
  # Source-data qualification supplies a private directory. Installed/release
  # tests leave this unset and exercise only the exact package's bundled data.
  load.example <- function(name) {
    env <- new.env(parent = baseenv())
    source <- Sys.getenv("NP_TEST_BUNDLED_DATA_SOURCE", "")
    if (nzchar(source)) load(file.path(source, paste0(name, ".rda")), env)
    else data(list = name, package = pkg, envir = env)
    env
  }
  wage <- load.example("wage1")
  oecd <- load.example("oecdpanel")
  expected <- list(basis = "glp", degree = rep(1L, 4L),
    bernstein.basis = FALSE, regtype.engine = "lp", basis.engine = "glp",
    degree.engine = rep(1L, 4L), bernstein.basis.engine = FALSE)
  for (name in names(expected))
    expect_identical(oecd$bw[[name, exact = TRUE]], expected[[name]], info = name)
  # Replay the existing wage example's explicit training-subset recipe.
  # This does not assert that the original historical search retained row IDs.
  ii <- sample(seq_len(nrow(wage$wage1)), replace = FALSE)
  cases <- list(
    all = list(bws = wage$bw.all, data = wage$wage1),
    subset = list(bws = wage$bw.subset, data = wage$wage1[ii[seq_len(400L)], , drop = FALSE]),
    oecd = list(bws = oecd$bw, data = oecd$oecdpanel))
  b <- oecd$bw
  reference.bw <- npregbw(b$formula, data = oecd$oecdpanel, bws = b$bw,
    bandwidth.compute = FALSE, regtype = b$regtype, bwmethod = b$method,
    bwtype = b$type, ckertype = b$ckertype, ckerorder = b$ckerorder,
    ukertype = b$ukertype, okertype = b$okertype)
  reference <- npreg(reference.bw, data = oecd$oecdpanel, se = TRUE)
  local_mocked_bindings(npregbw = function(...) stop("unexpected bandwidth selection"),
                       .package = pkg)
  for (name in names(cases)) {
    entry <- cases[[name]]; before <- serialize(entry$bws, NULL); rng <- .Random.seed
    fit <- npreg(entry$bws, data = entry$data, se = TRUE)
    expect_true(all(is.finite(fitted(fit))), info = name)
    expect_true(all(is.finite(se(fit))), info = name)
    expect_identical(fit$bws$bw, entry$bws$bw, info = name)
    prediction <- predict(fit, newdata = entry$data[seq_len(5L), , drop = FALSE])
    expect_length(prediction, 5L)
    expect_true(all(is.finite(prediction)), info = name)
    plotted <- plot(fit, output = "data", errors = "none", neval = 5L)
    expect_s3_class(plotted, "np_plot_data")
    expect_gt(length(plotted$panels), 0L)
    for (panel in seq_along(plotted$panels)) {
      flat <- as.data.frame(plotted, panel = panel)
      expect_gt(nrow(flat), 0L)
      expect_true(is.numeric(flat$estimate))
      expect_true(all(is.finite(flat$estimate)), info = name)
    }
    frame <- model.frame(entry$bws$formula, data = entry$data)
    expect_length(plotted$panels, ncol(frame) - 1L)
    direct <- plot(entry$bws, xdat = frame[-1L], ydat = model.response(frame),
                   output = "data", errors = "none", neval = 5L)
    expect_s3_class(direct, "np_plot_data")
    expect_length(direct$panels, length(plotted$panels))
    for (panel in seq_along(direct$panels)) {
      direct.data <- as.data.frame(direct, panel = panel)
      fitted.data <- as.data.frame(plotted, panel = panel)
      expect_identical(direct.data, fitted.data, info = name)
    }
    expect_identical(serialize(entry$bws, NULL), before, info = name)
    expect_identical(.Random.seed, rng, info = name)
    if (name == "oecd") {
      expect_identical(fitted(fit), fitted(reference))
      expect_identical(se(fit), se(reference))
      expect_identical(prediction, predict(reference, newdata = entry$data[seq_len(5L), , drop = FALSE]))
    }
  }
})
