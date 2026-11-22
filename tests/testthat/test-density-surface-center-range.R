test_that("density surface boxes enclose the endpoints actually drawn", {
  old <- options(np.messages = FALSE); on.exit(options(old), add = TRUE)
  set.seed(42); d <- data.frame(x = rnorm(90), z = rnorm(90))
  f <- npudens(tdat = d, bws = c(.45, .55))
  draw <- function(center, band, zlim = NULL, fit = f) {
    grDevices::pdf(NULL); on.exit(grDevices::dev.off(), add = TRUE)
    grDevices::dev.control(displaylist = "enable")
    args <- list(fit, errors = "asymptotic", center = center, view = "fixed",
                 neval = 7L, output = "plot-data", zlim = zlim)
    if (band != "default") args$band <- band
    if (band %in% c("simultaneous", "all")) args$gaussian.draws <- 199L
    seed <- .Random.seed
    out <- suppressWarnings(do.call(plot, args))
    expect_identical(.Random.seed, seed)
    display <- grDevices::recordPlot()[[1L]]
    boxes <- lapply(display, function(e) {
      z <- e[[2L]]
      if (identical(tryCatch(z[[1L]]$name, error = function(e) ""), "C_persp"))
        list(limits = as.double(z[[7L]]), background = all(is.na(z[[4L]]))) else NULL
    })
    boxes <- Filter(Negate(is.null), boxes)
    expect_identical(vapply(boxes, "[[", logical(1), "background"), c(TRUE, FALSE))
    expect_identical(boxes[[1L]]$limits, boxes[[2L]]$limits)
    list(out = out, box = boxes[[2L]]$limits)
  }
  for (center in c("estimate", "bias-corrected")) {
    for (band in c("default", "pointwise", "bonferroni", "simultaneous", "all")) {
      p <- draw(center, band)
      panel <- p$out$panels[[1L]]
      expected <- range(c(panel$values$estimate, panel$values$bias_corrected,
                          panel$intervals$lower, panel$intervals$upper), finite = TRUE)
      expect_true(p$box[1L] <= expected[1L] + 1e-12, info = paste(center, band))
      expect_true(p$box[2L] >= expected[2L] - 1e-12, info = paste(center, band))
    }
  }
  # Explicit limits remain a user decision even if they intentionally clip.
  explicit <- c(0, .1)
  expect_identical(draw("bias-corrected", "pointwise", explicit)$box, explicit)
  x <- data.frame(x = c(seq(.01,.19,length.out=25),seq(.81,.99,length.out=25)))
  dd <- data.frame(x=x$x, z=rev(x$x))
  sparse <- npudens(tdat=dd, bws=c(.08,.08), ckertype="uniform")
  for (band in c("pointwise", "simultaneous")) {
    p<-draw("estimate", band, fit=sparse)
    expect_true(all(is.finite(p$box)))
    expect_true(any(p$out$panels[[1L]]$intervals$status != "available"))
  }
})
