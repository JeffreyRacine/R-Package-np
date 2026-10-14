# Small, fixed-bandwidth contracts for plotting-only repairs. The affine test
# supplies known asymmetric intervals so it does not depend on bootstrap noise.
plot_audit_method <- function(object) {
  ns <- asNamespace(getNamespaceName(environment(npudens)))
  method <- getS3method("plot", class(object)[1L], envir = ns)
  stopifnot(identical(environment(method), ns))
  method
}

plot_audit_device <- function() {
  file <- tempfile(fileext = ".pdf")
  grDevices::pdf(file)
  device <- grDevices::dev.cur()
  withr::defer({
    if (!is.null(grDevices::dev.list()) && device %in% grDevices::dev.list())
      grDevices::dev.off(device)
    unlink(file)
  }, envir = parent.frame())
}

test_that("quantile plots reject unsupported derivative orders without narrowing LSQ", {
  withr::local_options(np.messages = FALSE, np.plot.progress = FALSE)
  withr::local_seed(137L)
  x <- data.frame(x = seq(-1, 1, length.out = 36L))
  y <- x$x^2 + rnorm(nrow(x), sd = .25)
  b <- npcdistbw(xdat = x, ydat = y, bws = c(.3, .5),
    bandwidth.compute = FALSE, regtype = "lp", degree = 2L)
  q <- npqreg(b, txdat = x, tydat = y, tau = .5)
  p <- plot_audit_method(q)
  first <- p(q, gradients = TRUE, gradient.order = 1L, errors = "none",
    output = "data", neval = 3L)
  expect_true(all(is.finite(first[[1L]]$quantgrad)))
  expect_error(p(q, gradients = TRUE, gradient.order = 2L, errors = "none",
    output = "data", neval = 3L), "conditional-quantile plots support only first derivatives",
    fixed = TRUE)
  expect_error(plot_audit_method(b)(b, xdat = x, ydat = y, quantreg = TRUE,
    gradients = TRUE, gradient.order = 2L, errors = "none", output = "data",
    neval = 3L), "conditional-quantile plots support only first derivatives", fixed = TRUE)
  l <- nplsqreg(bws = .5, txdat = x, tydat = y, tau = .5,
    bandwidth.compute = FALSE, regtype = "lp", degree = 2L)
  lp <- plot_audit_method(l)
  first <- lp(l, gradients = TRUE, gradient.order = 1L, errors = "none",
    output = "data", neval = 3L)
  second <- lp(l, gradients = TRUE, gradient.order = 2L, errors = "none",
    output = "data", neval = 3L)
  expect_true(all(is.finite(second[[1L]]$quantgrad)))
  expect_false(identical(first[[1L]]$quantgrad, second[[1L]]$quantgrad))
})

test_that("data-only plots preserve graphics while renderer exits clean up", {
  withr::local_options(np.messages = FALSE, np.plot.progress = FALSE)
  withr::local_seed(139L)
  x <- data.frame(x = seq(-1, 1, length.out = 32L))
  y <- .4 * x$x + rnorm(nrow(x), sd = .3)
  b <- npcdistbw(xdat = x, ydat = y, bws = c(.3, .5), bandwidth.compute = FALSE)
  d <- npudens(npudensbw(dat = x, bws = .5, bandwidth.compute = FALSE), tdat = x)
  rb <- npregbw(xdat = x, ydat = y, bws = .5, bandwidth.compute = FALSE)
  r <- npreg(rb, txdat = x, tydat = y)
  fits <- list(d, r, npcdist(b, txdat = x, tydat = y), npqreg(b, txdat = x, tydat = y))
  plot_audit_device()
  reset <- function() {
    graphics::plot.new()
    graphics::plot.window(c(2, 7), c(-3, 4))
    graphics::par(new = TRUE)
  }
  for (fit in fits) {
    reset()
    before <- graphics::par(c("new", "mfrow", "mfg", "cex", "usr"))
    out <- plot_audit_method(fit)(fit, errors = "none", output = "data",
      perspective = FALSE, neval = 3L)
    expect_type(out, "list")
    expect_identical(graphics::par(names(before)), before)
  }
  # Fixed regression must also preserve caller-owned state when data validation
  # fails, both before and after the plot options have been normalized.
  graphics::par(mfrow = c(2, 2))
  for (expr in list(
    quote(plot_audit_method(rb)(rb, xdat = x, errors = "none", output = "data")),
    quote(plot_audit_method(r)(r, errors = "none", output = "data", xtrim = 2)))) {
    reset()
    before <- graphics::par(c("new", "mfrow", "mfg", "cex", "usr"))
    expect_error(eval(expr))
    expect_identical(graphics::par(names(before)), before)
  }
  for (fit in list(d, r)) {
    reset()
    before <- graphics::par(c("mfrow", "cex", "usr"))
    plot_audit_method(fit)(fit, errors = "none", output = "plot-data", neval = 3L)
    expect_identical(graphics::par("new"), FALSE)
    expect_identical(graphics::par(c("mfrow", "cex")), before[c("mfrow", "cex")])
    expect_false(identical(graphics::par("usr"), before$usr))
    reset()
    before <- graphics::par(c("mfrow", "cex"))
    expect_error(plot_audit_method(fit)(fit, errors = "none", output = "plot-data",
      neval = 3L, col = "not-a-valid-colour-contract-probe"), "invalid colou?r name")
    expect_identical(graphics::par("new"), FALSE)
    expect_identical(graphics::par(names(before)), before)
  }
})

test_that("conmode pmzsd gradient intervals use the pointwise renderer", {
  withr::local_options(np.messages = FALSE, np.plot.progress = FALSE)
  withr::local_seed(13841L)
  x <- data.frame(x = seq(-1.2, 1.2, length.out = 40L))
  y <- factor(ifelse(runif(nrow(x)) < plogis(1.2 * x$x), "north", "south"),
    levels = c("south", "north"))
  b <- npcdensbw(xdat = x, ydat = y, bws = c(.18, .52),
    bandwidth.compute = FALSE, regtype = "ll")
  fit <- npconmode(b, txdat = x, tydat = y, probabilities = TRUE,
    gradients = TRUE, level = "north", se = TRUE)
  p <- plot_audit_method(fit)
  plot_audit_device()
  for (view in c("fixed", "sample")) {
    pw <- p(fit, gradients = TRUE, level = "north", view = view, neval = 3L,
      errors = "asymptotic", band = "pointwise", output = "plot-data", layout = FALSE)
    pm <- p(fit, gradients = TRUE, level = "north", view = view, neval = 3L,
      errors = "asymptotic", band = "pmzsd", output = "plot-data", layout = FALSE)
    expect_true(all(is.finite(pw[[1L]]$lower)))
    expect_identical(pm, pw)
  }
})

test_that("negative index coefficients exchange asymmetric all-band endpoints", {
  withr::local_options(np.messages = FALSE, np.plot.progress = FALSE)
  withr::local_seed(27183L)
  n <- 32L
  x <- data.frame(x1 = runif(n, -1, 1), x2 = runif(n, -1, 1))
  beta <- c(1, -.8)
  y <- sin(1.4 * drop(as.matrix(x) %*% beta)) + rnorm(n, sd = .2)
  b <- npindexbw(xdat = x, ydat = y, bws = c(beta, .55),
    bandwidth.compute = FALSE, method = "ichimura", regtype = "ll")
  fit <- npindex(b, txdat = x, tydat = y)
  pkg <- getNamespaceName(environment(npudens))
  draw.all <- getFromNamespace("draw.all.error.types", pkg)
  panels <- draws <- list()
  local_mocked_bindings(
    compute.bootstrap.errors = function(idx.eval, ...) {
      m <- nrow(idx.eval)
      err <- cbind(lower = seq_len(m) / 100, upper = seq_len(m) / 40, center = 0)
      list(boot.err = err, boot.all.err = list(pointwise = err[, 1:2],
        simultaneous = 1.2 * err[, 1:2], bonferroni = 1.4 * err[, 1:2]))
    },
    draw.all.error.types = function(ex, center, all.err, ...) {
      panels[[length(panels) + 1L]] <<- list(ex = ex, center = center, bands = all.err)
      draw.all(ex = ex, center = center, all.err = all.err, ...)
    },
    draw.errors = function(ex, ely, ehy, ...) {
      draws[[length(draws) + 1L]] <<- list(ex = ex, lower = ely, upper = ehy)
    }, .package = pkg)
  plot_audit_device()
  plot_audit_method(fit)(fit, gradients = TRUE, errors = "bootstrap",
    bootstrap = "inid", band = "all", B = 39L, neval = 4L,
    output = "plot-data", data.overlay = FALSE)
  expect_length(panels, 2L)
  expect_length(draws, 6L)
  expect_identical(panels[[1L]]$ex, panels[[2L]]$ex)
  expect_equal(panels[[2L]]$center, beta[2L] * panels[[1L]]$center, tolerance = 1e-12)
  for (i in seq_along(c("pointwise", "simultaneous", "bonferroni"))) {
    err <- panels[[1L]]$bands[[i]]
    lower <- panels[[1L]]$center - err[, 1L]
    upper <- panels[[1L]]$center + err[, 2L]
    expect_true(any(err[, 1L] != err[, 2L]))
    expect_equal(draws[[i]]$lower, lower, tolerance = 1e-12)
    expect_equal(draws[[i]]$upper, upper, tolerance = 1e-12)
    expect_equal(draws[[i + 3L]]$lower, beta[2L] * upper, tolerance = 1e-12)
    expect_equal(draws[[i + 3L]]$upper, beta[2L] * lower, tolerance = 1e-12)
  }
})
