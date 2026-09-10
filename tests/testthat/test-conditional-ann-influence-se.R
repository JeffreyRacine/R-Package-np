# Small full-row oracle only: production uncertainty must not retain n-by-m
# weights. No search, package loading or MPI pool initialization belongs here.
local({
  n <- 31L
  i <- seq_len(n)
  x <- data.frame(x1 = seq(-1, 1, length.out = n),
                  x2 = sin(i * .71) + i/100)
  y <- data.frame(yu = factor(rep(c("u", "v"), length.out = n)),
                  yo = ordered(rep(c("lo", "mid", "hi"), length.out = n),
                               levels = c("lo", "mid", "hi")))
  ex <- x[c(4L, 9L, 21L, 28L), , drop = FALSE]
  ex$x1 <- ex$x1 + .031
  ey <- y[c(7L, 11L, 19L, 23L), , drop = FALSE]
  lu <- .25
  lo <- .35

  # Recover complete normalized donor rows through the independent public
  # kernel-sum interface, not the conditional estimator's SE implementation.
  rows <- function(h, op, kernel) {
    ans <- npksum(bws = h, txdat = x, exdat = ex, tydat = diag(n),
      bandwidth.divide = TRUE, operator = op, bwtype = "adaptive_nn",
      ckertype = kernel)$ksum
    expect_identical(dim(ans), c(nrow(ex), n))
    t(ans)
  }

  # The ordered Y operator is finite-support Racine-Li-Yan, normalized over
  # the fixture's three observed categories. For a categorical CDF, sum its
  # category probabilities through the query category.
  response <- function(cdf) {
    vapply(seq_len(nrow(ey)), function(j) {
      vapply(seq_len(n), function(donor) {
        p <- lo^abs(as.integer(y$yo[donor]) - seq_len(nlevels(y$yo)))
        p <- p/sum(p)
        if (cdf) {
          sum(p[seq_len(as.integer(ey$yo[j]))])
        } else {
          unordered <- if (y$yu[donor] == ey$yu[j]) 1-lu else
            lu/(nlevels(y$yu)-1L)
          unordered * p[as.integer(ey$yo[j])]
        }
      }, numeric(1))
    }, numeric(n))
  }

  cells <- data.frame(cdf = c(FALSE, FALSE, TRUE, TRUE),
    kernel = c("gaussian", "epanechnikov", "epanechnikov", "gaussian"),
    regtype = c("lc", "lp", "lc", "lp"),
    tree = c(FALSE, TRUE, TRUE, FALSE), extended = c(FALSE, FALSE, FALSE, TRUE))

  run <- function() {
    for (cell in seq_len(nrow(cells))) {
      spec <- cells[cell, ]
      label <- paste(if (spec$cdf) "CDF" else "density", spec$kernel,
                     if (spec$regtype == "lp") "LP-zero" else "LC",
                     if (spec$extended) "extended" else "ordinary")
      test_that(paste("conditional ANN influence SEs preserve", label), {
        old <- options(np.messages = FALSE, np.tree = spec$tree,
                       np.extendednn = TRUE)
        on.exit(options(old), add = TRUE)
        h <- if (spec$extended) rep(3L*n, 2L) else c(9L, 11L)
        yy <- if (spec$cdf) y["yo"] else y
        eyy <- if (spec$cdf) ey["yo"] else ey
        ybw <- if (spec$cdf) lo else c(lu, lo)
        bwfun <- if (spec$cdf) npcdistbw else npcdensbw
        fitfun <- if (spec$cdf) npcdist else npcdens
        bargs <- list(xdat = x, ydat = yy, bws = c(ybw, h),
          bandwidth.compute = FALSE, bwtype = "adaptive_nn",
          cxkertype = spec$kernel, oykertype = "racineliyan",
          regtype = spec$regtype)
        if (spec$regtype == "lp") bargs$degree <- c(0L, 0L)
        bw <- do.call(bwfun, bargs)
        fargs <- list(bws = bw, txdat = x, tydat = yy, exdat = ex,
                      eydat = eyy, gradients = TRUE)
        off <- do.call(fitfun, c(fargs, list(se = FALSE)))
        on <- do.call(fitfun, c(fargs, list(se = TRUE)))
        expect_identical(fitted(on), fitted(off))
        expect_identical(gradients(on), gradients(off))
        expect_null(off[["conderr", exact = TRUE]])
        expect_null(off[["congerr", exact = TRUE]])

        w <- rows(h, rep("normal", 2L), spec$kernel)
        z <- response(spec$cdf)
        den <- colSums(w)
        mean <- colSums(w*z)/den
        a <- sweep(w, 2L, den, "/")
        residual <- sweep(z, 2L, mean, "-")
        level.se <- sqrt(n/(n-1L) * colSums((a*residual)^2))
        expect_true(all(is.finite(den) & den > 0))
        expect_equal(as.numeric(fitted(on)), mean, tolerance = 5e-11)
        expect_equal(as.numeric(se(on)), level.se, tolerance = 5e-11)

        derivative <- derivative.se <- matrix(NA_real_, nrow(ex), 2L)
        for (j in seq_len(2L)) {
          op <- rep("normal", 2L)
          op[j] <- "derivative"
          dw <- rows(h, op, spec$kernel)
          dden <- colSums(dw)
          derivative[, j] <- colSums(dw*z)/den - mean*dden/den
          # Differentiate the complete a*(z-mean): both denominator and
          # conditional-mean derivatives matter, even with fixed categorical Y.
          da <- sweep(dw, 2L, den, "/") -
            sweep(a, 2L, dden/den, "*")
          influence <- da*residual - sweep(a, 2L, derivative[, j], "*")
          derivative.se[, j] <- sqrt(n/(n-1L) * colSums(influence^2))
        }
        expect_equal(unname(gradients(on)), derivative, tolerance = 5e-11)
        expect_equal(unname(gradients(on, se = TRUE)), derivative.se,
                     tolerance = 5e-11)

        fargs$gradients <- FALSE
        levels.only <- do.call(fitfun, c(fargs, list(se = TRUE)))
        expect_equal(as.numeric(fitted(levels.only)), mean, tolerance = 5e-11)
        expect_equal(as.numeric(se(levels.only)), level.se, tolerance = 5e-11)
        if (spec$regtype == "lp") {
          bargs$regtype <- "lc"
          bargs$degree <- NULL
          fargs$bws <- do.call(bwfun, bargs)
          fargs$gradients <- TRUE
          lc <- do.call(fitfun, c(fargs, list(se = TRUE)))
          expect_equal(fitted(on), fitted(lc), tolerance = 5e-11)
          expect_equal(gradients(on), gradients(lc), tolerance = 5e-11)
          expect_equal(se(on), se(lc), tolerance = 5e-11)
          expect_equal(gradients(on, se = TRUE), gradients(lc, se = TRUE),
                       tolerance = 5e-11)
        }
      })
    }
  }

  # Use the current estimator namespace, never whichever package happens to
  # be installed or loaded elsewhere. MPI standard checks use its existing
  # local owner and do not spawn a pool; pooled execution is qualified apart.
  package <- getNamespaceName(environment(npcdens))
  if (package == "npRmpi")
    getFromNamespace(".npRmpi_with_local_regression", package)(run())
  else
    run()
})
