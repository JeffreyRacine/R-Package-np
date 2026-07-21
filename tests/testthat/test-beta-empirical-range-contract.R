.beta_half_spacing_bounds <- function(x) {
  u <- sort(unique(as.double(x)))
  c(
    lower = u[[1L]] - (u[[2L]] - u[[1L]]) / 2,
    upper = u[[length(u)]] +
      (u[[length(u)]] - u[[length(u) - 1L]]) / 2
  )
}

test_that("empirical-range beta operators equal explicit half-spacing bounds", {
  training <- data.frame(x = c(-2, -1.86, -1.3, -0.35, 0.8, 2.2, 3))
  evaluation <- data.frame(x = seq(-2, 3, length.out = 11L))
  bounds <- .beta_half_spacing_bounds(training$x)

  for (order in c(2L, 4L, 6L, 8L)) {
    for (bwtype in c("fixed", "generalized_nn", "adaptive_nn")) {
      bandwidth <- if (identical(bwtype, "fixed")) 0.55 else 3
      for (operator in c("normal", "integral", "convolution")) {
        common <- list(
          txdat = training, exdat = evaluation, bws = bandwidth,
          bwtype = bwtype, ckertype = "beta", ckerorder = order,
          operator = operator, return.kernel.weights = TRUE
        )
        empirical <- do.call(npksum, c(common, list(ckerbound = "range")))
        explicit <- do.call(npksum, c(common, list(
          ckerbound = "fixed", ckerlb = bounds[["lower"]],
          ckerub = bounds[["upper"]]
        )))

        expect_equal(empirical$kw, explicit$kw, tolerance = 2e-12)
        expect_equal(empirical$ksum, explicit$ksum, tolerance = 2e-12)
      }
    }
  }
})

test_that("empirical-range beta metadata and unconditional estimators are exact", {
  training <- data.frame(x = c(-2, -1.86, -1.3, -0.35, 0.8, 2.2, 3))
  evaluation <- data.frame(x = c(-2, -1.1, 0.4, 1.7, 3))
  response <- c(0.2, 0.5, -0.1, 0.8, 1.4, 1.1, 2.3)
  bounds <- .beta_half_spacing_bounds(training$x)
  range.args <- list(
    bws = 0.55, bandwidth.compute = FALSE,
    ckertype = "beta", ckerorder = 6, ckerbound = "range"
  )
  fixed.args <- modifyList(range.args, list(
    ckerbound = "fixed", ckerlb = bounds[["lower"]],
    ckerub = bounds[["upper"]]
  ))

  density.range.bw <- do.call(npudensbw, c(list(dat = training), range.args))
  density.fixed.bw <- do.call(npudensbw, c(list(dat = training), fixed.args))
  distribution.range.bw <- do.call(npudistbw, c(list(dat = training), range.args))
  distribution.fixed.bw <- do.call(npudistbw, c(list(dat = training), fixed.args))
  regression.range.bw <- do.call(npregbw, c(list(
    xdat = training, ydat = response, regtype = "lc"
  ), range.args))
  regression.fixed.bw <- do.call(npregbw, c(list(
    xdat = training, ydat = response, regtype = "lc"
  ), fixed.args))

  expect_identical(density.range.bw$ckerbound, "range")
  expect_equal(density.range.bw$ckerlb[density.range.bw$icon],
               bounds[["lower"]])
  expect_equal(density.range.bw$ckerub[density.range.bw$icon],
               bounds[["upper"]])
  expect_equal(fitted(npudens(
    bws = density.range.bw, tdat = training, edat = evaluation
  )), fitted(npudens(
    bws = density.fixed.bw, tdat = training, edat = evaluation
  )), tolerance = 2e-12)
  expect_equal(fitted(npudist(
    bws = distribution.range.bw, tdat = training, edat = evaluation
  )), fitted(npudist(
    bws = distribution.fixed.bw, tdat = training, edat = evaluation
  )), tolerance = 2e-12)
  expect_equal(fitted(npreg(
    bws = regression.range.bw, txdat = training, tydat = response,
    exdat = evaluation
  )), fitted(npreg(
    bws = regression.fixed.bw, txdat = training, tydat = response,
    exdat = evaluation
  )), tolerance = 2e-12)

  padded.endpoint <- npudist(
    bws = distribution.range.bw, tdat = training,
    edat = data.frame(x = bounds)
  )
  raw.endpoint <- npudist(
    bws = distribution.range.bw, tdat = training,
    edat = data.frame(x = range(training$x))
  )
  expect_equal(fitted(padded.endpoint), c(0, 1), tolerance = 0)
  expect_true(fitted(raw.endpoint)[[1L]] > 0)
  expect_true(fitted(raw.endpoint)[[2L]] < 1)
})

test_that("empirical-range beta gradients are finite at raw extrema", {
  training <- data.frame(x = c(-2, -1.86, -1.3, -0.35, 0.8, 2.2, 3))
  evaluation <- data.frame(x = c(-2, -0.4, 3))
  response <- c(0.2, 0.5, -0.1, 0.8, 1.4, 1.1, 2.3)
  bounds <- .beta_half_spacing_bounds(training$x)

  for (order in c(2L, 4L, 6L, 8L)) {
    for (bwtype in c("fixed", "generalized_nn", "adaptive_nn")) {
      bandwidth <- if (identical(bwtype, "fixed")) 0.55 else 3
      common <- list(
        txdat = training, tydat = response, exdat = evaluation,
        bws = bandwidth, bwtype = bwtype, gradients = TRUE,
        regtype = "lc", ckertype = "beta", ckerorder = order
      )
      expect_warning(
        empirical <- do.call(
          npreg, c(common, list(ckerbound = "range"))
        ),
        NA
      )
      explicit <- do.call(npreg, c(common, list(
        ckerbound = "fixed", ckerlb = bounds[["lower"]],
        ckerub = bounds[["upper"]]
      )))

      expect_identical(fitted(empirical), fitted(explicit))
      expect_identical(gradients(empirical), gradients(explicit))
      expect_identical(empirical$gerr, explicit$gerr)
      expect_true(all(is.finite(gradients(empirical)[, 1L])))
      expect_true(all(is.finite(empirical$gerr[, 1L])))
    }
  }

  constant <- npreg(
    bws = 0.55, txdat = training, tydat = rep(7, nrow(training)),
    exdat = evaluation[c(1L, 3L), , drop = FALSE], gradients = TRUE,
    regtype = "lc", ckertype = "beta", ckerorder = 8,
    ckerbound = "range"
  )
  expect_identical(as.double(gradients(constant)), c(0, 0))
  expect_identical(as.double(constant$gerr), c(0, 0))
})

test_that("empirical-range conditional beta resolves X and Y independently", {
  training.x <- data.frame(x = c(-2, -1.7, -1.05, -0.2, 0.9, 2.1, 3))
  training.y <- data.frame(y = c(10, 10.3, 11.1, 12.4, 13.8, 15.2, 16))
  evaluation.x <- data.frame(x = c(-2, -1, 0.5, 2, 3))
  evaluation.y <- data.frame(y = c(10, 10.8, 12.9, 15.4, 16))
  xbounds <- .beta_half_spacing_bounds(training.x$x)
  ybounds <- .beta_half_spacing_bounds(training.y$y)

  for (bwtype in c("fixed", "generalized_nn", "adaptive_nn")) {
    xbw <- if (identical(bwtype, "fixed")) 0.52 else 3
    ybw <- if (identical(bwtype, "fixed")) 0.66 else 3
    range.common <- list(
      xdat = training.x, ydat = training.y, bws = c(ybw, xbw),
      bandwidth.compute = FALSE, bwtype = bwtype,
      cxkertype = "beta", cxkerorder = 4, cxkerbound = "range",
      cykertype = "beta", cykerorder = 8, cykerbound = "range"
    )
    fixed.common <- modifyList(range.common, list(
      cxkerbound = "fixed", cxkerlb = xbounds[["lower"]],
      cxkerub = xbounds[["upper"]],
      cykerbound = "fixed", cykerlb = ybounds[["lower"]],
      cykerub = ybounds[["upper"]]
    ))
    density.range.bw <- do.call(npcdensbw, range.common)
    density.fixed.bw <- do.call(npcdensbw, fixed.common)
    distribution.range.bw <- do.call(npcdistbw, range.common)
    distribution.fixed.bw <- do.call(npcdistbw, fixed.common)

    expect_identical(density.range.bw$cxkerbound, "range")
    expect_identical(density.range.bw$cykerbound, "range")
    expect_equal(density.range.bw$cxkerlb[density.range.bw$ixcon],
                 xbounds[["lower"]])
    expect_equal(density.range.bw$cxkerub[density.range.bw$ixcon],
                 xbounds[["upper"]])
    expect_equal(density.range.bw$cykerlb[density.range.bw$iycon],
                 ybounds[["lower"]])
    expect_equal(density.range.bw$cykerub[density.range.bw$iycon],
                 ybounds[["upper"]])

    density.range <- npcdens(
      bws = density.range.bw, txdat = training.x, tydat = training.y,
      exdat = evaluation.x, eydat = evaluation.y
    )
    density.fixed <- npcdens(
      bws = density.fixed.bw, txdat = training.x, tydat = training.y,
      exdat = evaluation.x, eydat = evaluation.y
    )
    distribution.range <- npcdist(
      bws = distribution.range.bw, txdat = training.x, tydat = training.y,
      exdat = evaluation.x, eydat = evaluation.y
    )
    distribution.fixed <- npcdist(
      bws = distribution.fixed.bw, txdat = training.x, tydat = training.y,
      exdat = evaluation.x, eydat = evaluation.y
    )

    expect_equal(fitted(density.range), fitted(density.fixed),
                 tolerance = 4e-10)
    expect_equal(se(density.range), se(density.fixed), tolerance = 4e-10)
    expect_equal(fitted(distribution.range), fitted(distribution.fixed),
                 tolerance = 4e-10)
    expect_equal(se(distribution.range), se(distribution.fixed),
                 tolerance = 4e-10)
  }
})

test_that("empirical-range bounds reach automatic beta objective searches", {
  x <- data.frame(x = c(-2, -1.7, -1.05, -0.2, 0.9, 2.1, 3))
  y <- data.frame(y = c(10, 10.3, 11.1, 12.4, 13.8, 15.2, 16))
  xbounds <- .beta_half_spacing_bounds(x$x)
  ybounds <- .beta_half_spacing_bounds(y$y)
  unconditional.common <- list(
    dat = x, ckertype = "beta", ckerorder = 2,
    bwmethod = "cv.ml", bwtype = "fixed", bwscaling = FALSE,
    nmulti = 1L, itmax = 6L, bwsolver = "powell",
    scale.factor.init = 0.5,
    scale.factor.init.lower = 0.2,
    scale.factor.init.upper = 1.2,
    scale.factor.search.lower = 0.1
  )
  unconditional.range <- do.call(npudensbw, c(
    unconditional.common, list(ckerbound = "range")
  ))
  unconditional.fixed <- do.call(npudensbw, c(
    unconditional.common,
    list(ckerbound = "fixed", ckerlb = xbounds[["lower"]],
         ckerub = xbounds[["upper"]])
  ))
  expect_equal(unconditional.range$bw, unconditional.fixed$bw,
               tolerance = 2e-12)
  expect_equal(unconditional.range$fval, unconditional.fixed$fval,
               tolerance = 2e-12)

  conditional.common <- list(
    xdat = x, ydat = y,
    cxkertype = "beta", cxkerorder = 2,
    cykertype = "beta", cykerorder = 2,
    bwmethod = "cv.ml", bwtype = "fixed", bwscaling = FALSE,
    nmulti = 1L, itmax = 6L, bwsolver = "powell",
    scale.factor.init = 0.5,
    scale.factor.init.lower = 0.2,
    scale.factor.init.upper = 1.5,
    scale.factor.search.lower = 0.1
  )
  conditional.range <- do.call(npcdensbw, c(
    conditional.common,
    list(cxkerbound = "range", cykerbound = "range")
  ))
  conditional.fixed <- do.call(npcdensbw, c(
    conditional.common,
    list(
      cxkerbound = "fixed", cxkerlb = xbounds[["lower"]],
      cxkerub = xbounds[["upper"]],
      cykerbound = "fixed", cykerlb = ybounds[["lower"]],
      cykerub = ybounds[["upper"]]
    )
  ))
  expect_equal(conditional.range$xbw, conditional.fixed$xbw,
               tolerance = 2e-12)
  expect_equal(conditional.range$ybw, conditional.fixed$ybw,
               tolerance = 2e-12)
  expect_equal(conditional.range$fval, conditional.fixed$fval,
               tolerance = 2e-12)
})

test_that("empirical-range beta preserves exact range failure contracts", {
  training <- data.frame(x = c(-2, -1, 0, 1, 3))

  expect_error(
    npudens(
      tdat = training, edat = data.frame(x = 4.01), bws = 0.5,
      ckertype = "beta", ckerbound = "range"
    ),
    "Evaluation data violate 'ckerbound' bounds",
    fixed = TRUE
  )
  expect_error(
    npksum(
      txdat = data.frame(x = rep(1, 4)), bws = 0.5,
      ckertype = "beta", ckerbound = "range"
    ),
    "half-spacing bounds require at least two distinct finite training values",
    fixed = TRUE
  )
  expect_error(
    npudensbw(
      dat = data.frame(x = c(-.Machine$double.xmax, 0,
                             .Machine$double.xmax)),
      bws = 1, bandwidth.compute = FALSE,
      ckertype = "beta", ckerbound = "range"
    ),
    "half-spacing bounds are not finite and strictly outside",
    fixed = TRUE
  )
})

test_that("beta half-spacing uses distinct extrema and is affine equivariant", {
  x <- c(-2, -2, -1.5, 0.25, 2, 3, 3, 3)
  bounds <- .beta_half_spacing_bounds(x)
  bw <- npudensbw(
    ~x, bws = 0.4, bandwidth.compute = FALSE,
    ckertype = "beta", ckerbound = "range"
  )
  expect_equal(bw$ckerlb[bw$icon], bounds[["lower"]])
  expect_equal(bw$ckerub[bw$icon], bounds[["upper"]])

  shifted <- 7 - 3 * x
  shifted.bounds <- .beta_half_spacing_bounds(shifted)
  shifted.bw <- npudensbw(
    ~shifted, bws = 1.2, bandwidth.compute = FALSE,
    ckertype = "beta", ckerbound = "range"
  )
  expect_equal(shifted.bw$ckerlb[shifted.bw$icon],
               shifted.bounds[["lower"]])
  expect_equal(shifted.bw$ckerub[shifted.bw$icon],
               shifted.bounds[["upper"]])
  expect_equal(shifted.bounds,
               c(lower = 7 - 3 * bounds[["upper"]],
                 upper = 7 - 3 * bounds[["lower"]]))
})

test_that("exact fixed and legacy stored beta bounds are not padded", {
  x <- data.frame(x = c(-2, -1.5, 0.25, 2, 3))
  fixed <- npudensbw(
    dat = x, bws = 0.4, bandwidth.compute = FALSE,
    ckertype = "beta", ckerbound = "fixed", ckerlb = -2, ckerub = 3
  )
  expect_identical(unname(fixed$ckerlb[fixed$icon]), -2)
  expect_identical(unname(fixed$ckerub[fixed$icon]), 3)

  legacy <- npRmpi:::untangle(x)
  legacy$all.min.next <- NULL
  legacy$all.max.prev <- NULL
  resolved <- npRmpi:::npKernelBoundsResolve(
    dati = legacy, varnames = "x", kerbound = "range",
    range.policy = "beta_half_spacing"
  )
  expect_identical(unname(resolved$lb[legacy$icon]), -2)
  expect_identical(unname(resolved$ub[legacy$icon]), 3)
})

test_that("non-beta empirical ranges remain exact", {
  x <- c(-2, -1.5, 0.25, 2, 3)
  for (kernel in c("gaussian", "epanechnikov", "uniform")) {
    bw <- npudensbw(
      ~x, bws = 0.4, bandwidth.compute = FALSE,
      ckertype = kernel, ckerbound = "range"
    )
    expect_identical(unname(bw$ckerlb[bw$icon]), min(x))
    expect_identical(unname(bw$ckerub[bw$icon]), max(x))
  }
})
