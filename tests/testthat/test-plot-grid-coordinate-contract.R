# Coordinates, rather than row counts alone, protect the plot grid contract.
test_that("density default and explicit counts reach the same continuous grid", {
  if (exists("spawn_mpi_slaves", mode = "function")) {
    if (!spawn_mpi_slaves()) skip("Could not initialize MPI context")
    on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  }
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  x <- data.frame(x = rep(0:4, 8L))
  fit <- npudens(tdat = x, bws = 0.6)
  saved <- serialize(fit, NULL)
  capture <- new.env(parent = emptyenv())
  original <- graphics::plot
  local_mocked_bindings(plot = function(x, y = NULL, ...) {
    if (is.numeric(x)) {
      capture$x <- as.numeric(x)
      capture$y <- as.numeric(y)
      original(x, y, ...)
    } else original(x, ...)
  }, .package = getNamespaceName(environment(npudens)))
  grDevices::pdf(file = tempfile(fileext = ".pdf"))
  on.exit(grDevices::dev.off(), add = TRUE)
  for (count in c(50L, 17L)) {
    args <- list(x = fit)
    if (count != 50L) args$neval <- count
    do.call(plot, args)
    expect_equal(capture$x, seq(0, 4, length.out = count))
    reference <- npudens(bws = fit$bws, tdat = x,
                         edat = data.frame(x = capture$x))
    expect_equal(capture$y, as.numeric(fitted(reference)), tolerance = 1e-12)
  }
  plot(fit, neval = 50L)
  expect_equal(capture$x, seq(0, 4, length.out = 50L))
  expect_identical(serialize(fit, NULL), saved)
})

test_that("conditional-mode numeric grids retain fractions regardless of storage", {
  if (exists("spawn_mpi_slaves", mode = "function")) {
    if (!spawn_mpi_slaves()) skip("Could not initialize MPI context")
    on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  }
  old <- options(np.messages = FALSE)
  on.exit(options(old), add = TRUE)
  x <- data.frame(x = rep(0:4, 8L), z = rep(0:3, each = 10L))
  y <- data.frame(y = factor(rep(c("no", "yes"), 20L)))
  fit <- npconmode(txdat = x, tydat = y, bws = c(.2, .6, .7),
                  probabilities = TRUE, gradients = TRUE, level = "yes", regtype = "lc")
  double <- x
  double[] <- lapply(double, as.double)
  reference <- npconmode(txdat = double, tydat = y, bws = fit$bws,
                        probabilities = TRUE, gradients = TRUE, level = "yes")
  for (count in c(50L, 17L)) {
    for (gradient in c(FALSE, TRUE)) {
      out <- plot(fit, view = "fixed", neval = count, output = "data",
                  level = "yes", gradients = gradient, xq = c(.5, .5))
      ref <- plot(reference, view = "fixed", neval = count, output = "data",
                  level = "yes", gradients = gradient, xq = c(.5, .5))
      expect_equal(out$x$x, seq(0, 4, length.out = count))
      expect_equal(out$z$x, seq(0, 3, length.out = count))
      expect_equal(out, ref, tolerance = 1e-12)
    }
  }
  surface <- plot(fit, perspective = TRUE, neval = 7L, output = "data", level = "yes")
  expect_equal(sort(unique(surface$surface$x1)), seq(0, 4, length.out = 7L))
  expect_equal(sort(unique(surface$surface$x2)), seq(0, 3, length.out = 7L))
  expect_equal(surface, plot(reference, perspective = TRUE, neval = 7L,
                            output = "data", level = "yes"), tolerance = 1e-12)
})

test_that("conditional-mode coordinate casting preserves categorical semantics", {
  cast <- getFromNamespace(".np_plot_conmode_cast_like", getNamespaceName(environment(npconmode)))
  base <- getFromNamespace(".np_plot_conmode_base_row", getNamespaceName(environment(npconmode)))
  expect_identical(cast(c(.25, 1.75), 0:2), c(.25, 1.75))
  for (ordered in c(FALSE, TRUE)) {
    f <- factor(c("b", "a"), levels = c("a", "b"), ordered = ordered)
    expect_identical(cast(c("a", "b"), f), f[2:1])
  }
  expect_equal(base(data.frame(x = 0:3), .5)$x, 1.5)
})

test_that("categorical density bootstrap propagates admitted-owner failures", {
  if (exists("spawn_mpi_slaves", mode = "function")) {
    if (!spawn_mpi_slaves()) skip("Could not initialize MPI context")
    on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  }
  pkg <- getNamespaceName(environment(npudens))
  old <- options(np.messages = FALSE, np.categorical.compress = TRUE)
  on.exit(options(old), add = TRUE)
  x <- data.frame(x = factor(rep(c("a", "b"), 20L)))
  bw <- npudensbw(dat = x, bws = .2, bandwidth.compute = FALSE)
  boot <- getFromNamespace(".np_inid_boot_from_ksum_unconditional", pkg)
  calls <- new.env(parent = emptyenv())
  calls$dense <- 0L
  local_mocked_bindings(
    .np_density_cat_profile_kernel_matrix = function(...) stop("profile witness"),
    .np_ksum_unconditional_operator_fixed = function(...) {
      calls$dense <- calls$dense + 1L
      stop("alternate owner must not execute")
    }, .package = pkg)
  expect_error(boot(x, x[1:2, , drop = FALSE], bw, 2L, "normal",
                    counts = matrix(1, 40L, 2L)), "profile witness", fixed = TRUE)
  expect_identical(calls$dense, 0L)
})

test_that("categorical density bootstrap retains explicit compression on and off", {
  if (exists("spawn_mpi_slaves", mode = "function")) {
    if (!spawn_mpi_slaves()) skip("Could not initialize MPI context")
    on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  }
  old <- options(np.messages = FALSE, np.categorical.compress = TRUE)
  on.exit(options(old), add = TRUE)
  boot <- getFromNamespace(".np_inid_boot_from_ksum_unconditional",
                           getNamespaceName(environment(npudens)))
  counts <- cbind(rep(1L, 40L), rep(c(0L, 2L), 20L))
  for (ordered in c(FALSE, TRUE)) {
    x <- data.frame(x = factor(rep(c("a", "b"), 20L), ordered = ordered))
    bw <- npudensbw(dat = x, bws = .2, bandwidth.compute = FALSE)
    options(np.categorical.compress = TRUE)
    profile <- boot(x, x[1:2, , drop = FALSE], bw, 2L, "normal", counts = counts)
    options(np.categorical.compress = FALSE)
    dense <- boot(x, x[1:2, , drop = FALSE], bw, 2L, "normal", counts = counts)
    expect_equal(profile, dense, tolerance = 1e-12)
  }
})
