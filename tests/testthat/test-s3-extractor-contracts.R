test_that("core fitted and se extract stored quantities without recomputation", {
  dens <- structure(list(dens = c(1, 2), derr = c(0.1, 0.2)), class = "npdensity")
  expect_equal(fitted(dens), dens$dens)
  expect_equal(se(dens), dens$derr)

  dist <- structure(list(dist = c(0.2, 0.8), derr = c(0.03, 0.04)), class = "npdistribution")
  expect_equal(fitted(dist), dist$dist)
  expect_equal(se(dist), dist$derr)

  reg <- structure(
    list(mean = c(2, 3), merr = c(0.2, 0.3), residuals = TRUE, resid = c(-0.1, 0.1)),
    class = "npregression"
  )
  expect_equal(fitted(reg), reg$mean)
  expect_equal(se(reg), reg$merr)
  expect_equal(residuals(reg), reg$resid)

  qreg <- structure(
    list(quantile = matrix(1:4, nrow = 2), quanterr = matrix(0.1, nrow = 2, ncol = 2)),
    class = "qregression"
  )
  expect_equal(fitted(qreg), qreg$quantile)
  expect_equal(quantile(qreg), qreg$quantile)
  expect_equal(se(qreg), qreg$quanterr)
})

test_that("conditional and semiparametric extractors preserve stored shapes", {
  cdens <- structure(
    list(condens = c(1, 2), conderr = c(0.1, 0.2), proper.applied = FALSE),
    class = "condensity"
  )
  expect_equal(fitted(cdens), cdens$condens)
  expect_equal(se(cdens), cdens$conderr)

  cdist <- structure(
    list(condist = c(0.1, 0.9), conderr = c(0.02, 0.03), proper.applied = FALSE),
    class = "condistribution"
  )
  expect_equal(fitted(cdist), cdist$condist)
  expect_equal(se(cdist), cdist$conderr)

  pl <- structure(
    list(mean = c(1, 3), merr = c(0.3, 0.4), residuals = TRUE, resid = c(0.1, -0.1)),
    class = "plregression"
  )
  expect_equal(fitted(pl), pl$mean)
  expect_equal(residuals(pl), pl$resid)

  si <- structure(
    list(
      mean = c(1, 2),
      merr = c(0.2, 0.3),
      residuals = TRUE,
      resid = c(0.1, 0.2),
      beta = c(1, -0.5),
      xnames = c("x1", "x2"),
      betavcov = diag(2)
    ),
    class = "singleindex"
  )
  expect_equal(fitted(si), si$mean)
  expect_equal(se(si), si$merr)
  expect_equal(residuals(si), si$resid)
  expect_equal(unname(coef(si)), si$beta)
  expect_equal(names(coef(si)), si$xnames)
  expect_equal(vcov(si), si$betavcov)

  sc <- structure(
    list(mean = c(2, 4), merr = c(0.1, 0.2), residuals = TRUE, resid = c(0.2, 0.3)),
    class = "smoothcoefficient"
  )
  expect_equal(fitted(sc), sc$mean)
  expect_equal(se(sc), sc$merr)
  expect_equal(residuals(sc), sc$resid)
})

test_that("extractors fail clearly when optional quantities were not stored", {
  qreg <- structure(list(quantgrad = NA, quantgerr = NA), class = "qregression")
  expect_error(gradients(qreg), "gradients are not available")
  expect_error(gradients(qreg, errors = TRUE), "gradient standard errors are not available")
  expect_error(gradients(qreg, errors = "yes"), "'errors' must be TRUE or FALSE", fixed = TRUE)

  conmode <- structure(list(), class = "conmode")
  expect_error(gradients(conmode), "class-probability gradients/effects are not available")
  expect_error(gradients(conmode, errors = TRUE), "gradient standard errors are not available")
  expect_error(gradients(conmode, errors = "yes"), "'errors' must be TRUE or FALSE", fixed = TRUE)

  cdens <- structure(
    list(conderr = c(0.1, 0.2), proper.applied = TRUE),
    class = "condensity"
  )
  expect_error(se(cdens), "standard errors are unavailable for repaired conditional densities")

  cdist <- structure(
    list(conderr = c(0.1, 0.2), proper.applied = TRUE),
    class = "condistribution"
  )
  expect_error(se(cdist), "standard errors are unavailable for repaired conditional distributions")
})

test_that("conmode modal-class extraction uses predict rather than a fake mode method", {
  ns <- asNamespace(getNamespaceName(environment(npconmode)))
  s3 <- getNamespaceInfo(ns, "S3methods")

  expect_false(exists("mode.conmode", envir = ns, inherits = FALSE))
  expect_false(any(s3[, 1L] == "mode" & s3[, 2L] == "conmode"))

  fit <- structure(list(conmode = factor(c("a", "b", "a"))), class = "conmode")
  expect_equal(predict(fit), fit$conmode)
  expect_identical(mode(fit), "list")
})
