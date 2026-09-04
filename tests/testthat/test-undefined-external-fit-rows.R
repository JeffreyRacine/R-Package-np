test_that("undefined external smooth-coefficient rows are NA, not penalties", {

  d <- data.frame(y = sin(1:24), x = (1:24)/24, z = (1:24)/24)
  for (reg in c("lc", "ll", "lp")) {
    b <- suppressWarnings(npscoefbw(y ~ x | z, data = d, bws = .06,
      bandwidth.compute = FALSE, regtype = reg,
      degree = if (reg == "lp") 2L else NULL, ckertype = "uniform"))
    warnings <- character()
    fit <- withCallingHandlers(
      npscoef(b, exdat = data.frame(x = c(.5, .5)),
        ezdat = data.frame(z = c(.5, 100)), se = TRUE, betas = TRUE),
      warning = function(e) {
        if (grepl("local fit is undefined", conditionMessage(e), fixed = TRUE))
          warnings <<- c(warnings, conditionMessage(e))
        else
          expect_match(conditionMessage(e), "ignoring kernel order", fixed = TRUE)
        invokeRestart("muffleWarning")
      })
    expect_length(warnings, 1L)
    expect_match(warnings, "1 row(s) (2)", fixed = TRUE)
    expect_match(warnings, "returning NA", fixed = TRUE)
    for (field in c("mean", "merr", "beta", "grad", "gerr")) {
      value <- fit[[field]]
      expect_true(all(is.finite(if (is.matrix(value)) value[1L, ] else value[1L])))
      expect_true(all(is.na(if (is.matrix(value)) value[2L, ] else value[2L])))
    }
    valid <- suppressWarnings(npscoef(b, exdat = data.frame(x = .5),
      ezdat = data.frame(z = .5), se = TRUE, betas = TRUE))
    for (field in c("mean", "merr", "beta", "grad", "gerr")) {
      value <- fit[[field]]
      expect_equal(as.double(if (is.matrix(value)) value[1L, ] else value[1L]),
                   as.double(valid[[field]]), tolerance = 1e-13)
    }
    base <- suppressWarnings(npscoef(b))
    predicted <- suppressWarnings(predict(base, newdata = data.frame(x = c(.5, .5),
                                                                     z = c(.5, 100))))
    expect_equal(as.double(predicted), fit$mean, tolerance = 1e-13)
  }
})

test_that("undefined required training and bootstrap fits still fail", {

  d <- data.frame(y = sin(1:24), x = (1:24)/24, z = (1:24)/24)
  b <- suppressWarnings(npscoefbw(y ~ x | z, data = d, bws = .01,
    bandwidth.compute = FALSE, regtype = "lc", ckertype = "uniform"))
  expect_error(suppressWarnings(npscoef(b, leave.one.out = TRUE)),
               "local fit is undefined", fixed = TRUE)
  exact <- getFromNamespace(".np_inid_boot_from_scoef_exact", "np")
  expect_error(suppressWarnings(exact(
    txdat = d["x"], ydat = d$y, tzdat = d["z"],
    exdat = data.frame(x = .5), ezdat = data.frame(z = 100),
    bws = b, B = 2L)), "local fit is undefined", fixed = TRUE)
  ## Failure is not detected from coefficient magnitude.
  huge <- d
  huge$y <- huge$y * 1e155
  fit <- suppressWarnings(npscoef(b, txdat = huge["x"], tydat = huge$y,
    tzdat = huge["z"], exdat = data.frame(x = .5), ezdat = data.frame(z = .5)))
  expect_true(all(is.finite(fit$mean)))
  expect_gt(abs(fit$mean), sqrt(.Machine$double.xmax))
})
