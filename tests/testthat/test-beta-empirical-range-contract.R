test_that("empirical-range beta operators equal explicit sample bounds", {
  training <- data.frame(x = c(-2, -1.86, -1.3, -0.35, 0.8, 2.2, 3))
  evaluation <- data.frame(x = seq(-2, 3, length.out = 11L))

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
          ckerbound = "fixed", ckerlb = -2, ckerub = 3
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
  range.args <- list(
    bws = 0.55, bandwidth.compute = FALSE,
    ckertype = "beta", ckerorder = 6, ckerbound = "range"
  )
  fixed.args <- modifyList(range.args, list(
    ckerbound = "fixed", ckerlb = -2, ckerub = 3
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
  expect_equal(density.range.bw$ckerlb[density.range.bw$icon], -2)
  expect_equal(density.range.bw$ckerub[density.range.bw$icon], 3)
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

  endpoint <- npudist(
    bws = distribution.range.bw, tdat = training,
    edat = data.frame(x = c(-2, 3))
  )
  expect_equal(fitted(endpoint), c(0, 1), tolerance = 0)
})

test_that("empirical-range beta gradients retain sample-boundary jumps", {
  training <- data.frame(x = c(-2, -1.86, -1.3, -0.35, 0.8, 2.2, 3))
  evaluation <- data.frame(x = c(-2, -0.4, 3))
  response <- c(0.2, 0.5, -0.1, 0.8, 1.4, 1.1, 2.3)

  for (order in c(2L, 4L, 6L, 8L)) {
    for (bwtype in c("fixed", "generalized_nn", "adaptive_nn")) {
      bandwidth <- if (identical(bwtype, "fixed")) 0.55 else 3
      common <- list(
        txdat = training, tydat = response, exdat = evaluation,
        bws = bandwidth, bwtype = bwtype, gradients = TRUE,
        regtype = "lc", ckertype = "beta", ckerorder = order
      )
      empirical.warnings <- capture_warnings(
        empirical <- do.call(
          npreg, c(common, list(ckerbound = "range"))
        )
      )
      expect_identical(
        empirical.warnings,
        paste0(
          "beta regression gradient produced 2 infinite endpoint value(s) ",
          "and 0 undefined cancellation(s)"
        )
      )
      explicit <- suppressWarnings(do.call(npreg, c(common, list(
        ckerbound = "fixed", ckerlb = -2, ckerub = 3
      ))))

      expect_identical(fitted(empirical), fitted(explicit))
      expect_identical(gradients(empirical), gradients(explicit))
      expect_identical(empirical$gerr, explicit$gerr)
      expect_true(all(is.infinite(gradients(empirical)[c(1L, 3L), 1L])))
      expect_true(is.finite(gradients(empirical)[2L, 1L]))
      expect_true(all(is.na(empirical$gerr[c(1L, 3L), 1L])))
      expect_true(is.finite(empirical$gerr[2L, 1L]))
    }
  }

  constant <- suppressWarnings(npreg(
    bws = 0.55, txdat = training, tydat = rep(7, nrow(training)),
    exdat = evaluation[c(1L, 3L), , drop = FALSE], gradients = TRUE,
    regtype = "lc", ckertype = "beta", ckerorder = 8,
    ckerbound = "range"
  ))
  expect_identical(as.double(gradients(constant)), c(0, 0))
  expect_identical(as.double(constant$gerr), c(0, 0))
})

test_that("empirical-range conditional beta resolves X and Y independently", {
  training.x <- data.frame(x = c(-2, -1.7, -1.05, -0.2, 0.9, 2.1, 3))
  training.y <- data.frame(y = c(10, 10.3, 11.1, 12.4, 13.8, 15.2, 16))
  evaluation.x <- data.frame(x = c(-2, -1, 0.5, 2, 3))
  evaluation.y <- data.frame(y = c(10, 10.8, 12.9, 15.4, 16))

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
      cxkerbound = "fixed", cxkerlb = -2, cxkerub = 3,
      cykerbound = "fixed", cykerlb = 10, cykerub = 16
    ))
    density.range.bw <- do.call(npcdensbw, range.common)
    density.fixed.bw <- do.call(npcdensbw, fixed.common)
    distribution.range.bw <- do.call(npcdistbw, range.common)
    distribution.fixed.bw <- do.call(npcdistbw, fixed.common)

    expect_identical(density.range.bw$cxkerbound, "range")
    expect_identical(density.range.bw$cykerbound, "range")
    expect_equal(density.range.bw$cxkerlb[density.range.bw$ixcon], -2)
    expect_equal(density.range.bw$cxkerub[density.range.bw$ixcon], 3)
    expect_equal(density.range.bw$cykerlb[density.range.bw$iycon], 10)
    expect_equal(density.range.bw$cykerub[density.range.bw$iycon], 16)

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
    list(ckerbound = "fixed", ckerlb = -2, ckerub = 3)
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
      cxkerbound = "fixed", cxkerlb = -2, cxkerub = 3,
      cykerbound = "fixed", cykerlb = 10, cykerub = 16
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
      tdat = training, edat = data.frame(x = 3.01), bws = 0.5,
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
    "Invalid bounds for 'ckerbound'",
    fixed = TRUE
  )
})
