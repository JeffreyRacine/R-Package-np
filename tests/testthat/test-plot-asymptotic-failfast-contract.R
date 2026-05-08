test_that("plot contract: single-index asymptotic consumer payloads are supported", {
  skip_if_not_installed("np")

  set.seed(921)
  n <- 60
  x1 <- runif(n)
  x2 <- runif(n)
  y.cont <- sin(2 * pi * x1) + 0.35 * x2 + rnorm(n, sd = 0.08)
  y.bin <- as.numeric(x1 + x2 + rnorm(n, sd = 0.2) > 1)
  xdat <- data.frame(x1 = x1, x2 = x2)

  reg.cfgs <- list(
    list(regtype = "lc"),
    list(regtype = "ll"),
    list(regtype = "lp", basis = "tensor", degree = 2L)
  )
  neval <- 50L

  for (method in c("ichimura", "kleinspady")) {
    y.use <- if (identical(method, "kleinspady")) y.bin else y.cont
    for (bwtype in c("fixed", "generalized_nn", "adaptive_nn")) {
      h <- if (identical(bwtype, "fixed")) 0.30 else 6L
      for (cfg in reg.cfgs) {
        bw.args <- list(
          xdat = xdat,
          ydat = y.use,
          bws = c(1, h, h),
          bandwidth.compute = FALSE,
          method = method,
          bwtype = bwtype,
          regtype = cfg$regtype
        )
        if (!is.null(cfg$basis)) {
          bw.args$basis <- cfg$basis
          bw.args$degree <- cfg$degree
        }

        bw <- do.call(npindexbw, bw.args)
        out <- suppressWarnings(plot(
          bw,
          xdat = xdat,
          ydat = y.use,
          neval = neval,
          output = "data",
          errors = "asymptotic",
          band = "pointwise",
          perspective = FALSE
        ))[[1]]

        expect_s3_class(out, "singleindex")
        expect_equal(length(out$mean), neval, info = paste(method, bwtype, cfg$regtype, "mean length"))
        expect_equal(dim(out$merr), c(neval, 2L), info = paste(method, bwtype, cfg$regtype, "merr shape"))
        expect_true(all(is.finite(out$mean)), info = paste(method, bwtype, cfg$regtype, "mean finite"))
        expect_true(all(is.finite(out$merr)), info = paste(method, bwtype, cfg$regtype, "merr finite"))
      }
    }
  }
})

test_that("plot contract: partially linear asymptotic consumer payloads are supported", {
  skip_if_not_installed("np")

  set.seed(922)
  n <- 55
  x <- runif(n)
  z <- runif(n)
  y <- sin(2 * pi * z) + 0.5 * x + rnorm(n, sd = 0.08)

  reg.cfgs <- list(
    list(regtype = "lc"),
    list(regtype = "ll"),
    list(regtype = "lp", basis = "tensor", degree = 2L)
  )

  for (bwmethod in c("cv.ls", "cv.aic")) {
    for (bwtype in c("fixed", "generalized_nn", "adaptive_nn")) {
      h <- if (identical(bwtype, "fixed")) 0.30 else 6L
      for (cfg in reg.cfgs) {
        bw.args <- list(
          xdat = data.frame(x = x),
          zdat = data.frame(z = z),
          ydat = y,
          bws = matrix(c(h, h), nrow = 2L),
          bandwidth.compute = FALSE,
          bwmethod = bwmethod,
          bwtype = bwtype,
          regtype = cfg$regtype
        )
        if (!is.null(cfg$basis)) {
          bw.args$basis <- cfg$basis
          bw.args$degree <- cfg$degree
        }

        bw <- do.call(npplregbw, bw.args)
        out <- suppressWarnings(plot(
          bw,
          xdat = data.frame(x = x),
          ydat = y,
          zdat = data.frame(z = z),
          output = "data",
          errors = "asymptotic",
          band = "pointwise",
          perspective = FALSE
        ))[[1]]

        expect_s3_class(out, "plregression")
        expect_true(length(out$mean) > 0L, info = paste(bwmethod, bwtype, cfg$regtype, "mean length"))
        expect_equal(dim(out$merr), c(length(out$mean), 2L), info = paste(bwmethod, bwtype, cfg$regtype, "merr shape"))
        expect_true(all(is.finite(out$mean)), info = paste(bwmethod, bwtype, cfg$regtype, "mean finite"))
        expect_true(all(is.finite(out$merr)), info = paste(bwmethod, bwtype, cfg$regtype, "merr finite"))
      }
    }
  }
})
