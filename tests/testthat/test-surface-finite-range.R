test_that("automatic all-band boxes retain finite estimates without intervals", {
  fun <- getFromNamespace(".np_plot_all_surface_range", "npRmpi")
  expect_identical(fun(c(-2, 5, NA), c(NA, 1), c(NA, 4), list(c(-1, NA)), list(c(3, NA))), c(-2, 5))
  expect_identical(fun(c(NA, NaN), NA_real_, Inf, list(-Inf), list(NA_real_)), c(-1, 1))
  expect_identical(fun(c(1, 2), c(0, 1), c(2, 3), list(-1), list(4)), c(-1, 4))
})

test_that("automatic surface boxes tolerate unavailable targets and preserve data", {
  skip_on_cran() # Actual renderer and family/kernel/band matrix.
  old <- options(np.messages = FALSE); on.exit(options(old), add = TRUE)
  set.seed(9782)
  d <- data.frame(x = c(runif(75, 0, .3), runif(75, .7, 1)))
  d$y <- d$x + rexp(150)
  draw <- function(fit, errors, band, zlim = NULL, neval = 12L) {
    grDevices::pdf(NULL); on.exit(grDevices::dev.off(), add = TRUE)
    grDevices::dev.control(displaylist = "enable")
    args <- list(fit, errors = errors, view = "fixed",
                 neval = neval, output = "plot-data", zlim = zlim)
    if (errors != "none") args$band <- band
    if (band == "all" && errors != "none") args$gaussian_draws <- 99L
    seed <- .Random.seed
    out <- suppressWarnings(do.call(plot, args))
    expect_identical(.Random.seed, seed)
    boxes <- Filter(Negate(is.null), lapply(grDevices::recordPlot()[[1L]], function(e) {
      z <- e[[2L]]
      if (identical(tryCatch(z[[1L]]$name, error = function(e) ""), "C_persp") &&
          any(is.finite(z[[4L]])))
        as.double(z[[7L]]) else NULL
    }))
    expect_length(boxes, 1L)
    panel <- out$panels[[1L]]
    if (is.null(zlim)) {
      expected <- range(c(panel$values$estimate, panel$intervals$lower,
                          panel$intervals$upper), finite = TRUE)
      expect_true(boxes[[1L]][1L] <= expected[1L] + 1e-12)
      expect_true(boxes[[1L]][2L] >= expected[2L] - 1e-12)
    } else expect_identical(boxes[[1L]], zlim)
    args$output <- "data"
    data <- suppressWarnings(do.call(plot, args))
    expect_identical(out$panels, data$panels)
    expect_identical(.Random.seed, seed)
    out
  }
  for (kernel in c("uniform", "epanechnikov")) {
    for (family in c("npcdens", "npcdist")) {
      bw <- do.call(get(paste0(family, "bw")), list(xdat = d["x"], ydat = d["y"],
        bws = c(.3, if (kernel == "uniform") .08 else .035), bandwidth.compute = FALSE, cxkertype = kernel, cykertype = kernel))
      fit <- do.call(get(family), list(bws = bw, txdat = d["x"], tydat = d["y"]))
      for (errors in c("none", "asymptotic")) {
        p <- draw(fit, errors, "pointwise")
        expect_true(anyNA(p$panels[[1L]]$values$estimate))
      }
      draw(fit, "asymptotic", "all")
      draw(fit, "none", "pointwise", c(-.1, 1.5))
    }
  }
  # A higher-order density can have a finite estimate where SE is unavailable.
  set.seed(3); g <- rbinom(200, 1, .5)
  dd <- data.frame(a = ifelse(g == 1, rnorm(200, -2, .3), rnorm(200, 2, .3)),
                   b = ifelse(g == 1, rnorm(200, 0, .3), rnorm(200, 3, .3)))
  fit <- npudens(~a + b, data = dd, bws = c(.6, .6), bandwidth.compute = FALSE,
                 ckertype = "epanechnikov", ckerorder = 4)
  draw(fit, "asymptotic", "all", neval = 11L)
})
