iv_interface_fixture <- function(n = 24L) {
  set.seed(9010)
  w <- rnorm(n)
  v <- rnorm(n, sd = 0.24)
  z <- 0.4 * w + v
  y <- z^2 - 0.45 * v + rnorm(n, sd = 0.06)
  data.frame(y = y, z = z, w = w)
}

iv_interface_package <- function() {
  if ("npRmpi" %in% loadedNamespaces()) "npRmpi" else "np"
}

test_that("IV formula grammar is explicit and transformation safe", {
  dat <- iv_interface_fixture()
  parse.iv <- getFromNamespace(".np_iv_parse_formula", iv_interface_package())
  parsed <- parse.iv(y ~ I(z + 0) | w, "npregiv()")
  expect_named(parsed$roles, c("z", "w"))
  expect_identical(parsed$roles$z$labels, "I(z + 0)")
  dat$wf <- ordered(ifelse(dat$w > 0, "high", "low"),
                    levels = c("low", "high"))
  parsed.factor <- parse.iv(y ~ z | wf, "npregiv()")
  frame.factor <- model.frame(parsed.factor$combined, data = dat)
  role.frame <- getFromNamespace(".np_iv_role_frame", iv_interface_package())
  expect_s3_class(role.frame(frame.factor, parsed.factor$roles$w,
                             "npregiv()")$wf, "ordered")

  expect_error(parse.iv(y ~ z:w | w, "npregiv()"),
               "interaction operators", fixed = TRUE)
  expect_error(parse.iv(y ~ . | w, "npregiv()"),
               "does not support '.'", fixed = TRUE)
  w <- dat$w
  expect_error(npregiv(y ~ z | w, data = dat[c("y", "z")],
                       nmulti = 1L, iterate.max = 2L),
               "data must contain columns", fixed = TRUE)
  z <- dat$z
  expect_error(npregiv(y ~ z | w, data = dat,
                       newdata = data.frame(other = seq_len(nrow(dat))),
                       nmulti = 1L, iterate.max = 2L),
               "newdata must contain columns", fixed = TRUE)
})

test_that("npregiv formula, modern degree spelling, and accessors preserve state", {
  if (identical(iv_interface_package(), "npRmpi"))
    skip("installed MPI route sentinel owns npRmpi numerical coverage")
  dat <- iv_interface_fixture()

  native <- suppressWarnings(npregiv(
    y = dat$y, z = dat$z, w = dat$w,
    nmulti = 1L, iterate.max = 2L
  ))
  formula <- suppressWarnings(npregiv(
    y ~ z | w, data = dat,
    nmulti = 1L, iterate.max = 2L
  ))
  modern.ll <- suppressWarnings(npregiv(
    y = dat$y, z = dat$z, w = dat$w, regtype = "ll",
    nmulti = 1L, iterate.max = 2L
  ))

  expect_identical(formula$phi, native$phi)
  expect_identical(modern.ll$phi, native$phi)
  expect_identical(formula$call[[1L]], quote(npregiv))
  expect_false("..." %in% names(formula$call))
  expect_null(attr(formula$call, ".Environment"))
  expect_identical(fitted(formula), formula$phi)
  expect_identical(gradients(formula), formula$phi.deriv.1)
  expect_equal(residuals(formula), dat$y - formula$phi, tolerance = 0)
  expect_s3_class(summary(formula), "summary.npregiv")
  expect_identical(summary(formula)$evaluated, length(formula$norm.stop))
  expect_identical(formula$smoothing.spec$effective$degree, 1L)
  expect_error(
    npregiv(y = dat$y, z = dat$z, w = dat$w,
            p = 1L, regtype = "lc", nmulti = 1L, iterate.max = 2L),
    "conflicting smoothing controls", fixed = TRUE
  )
  expect_error(
    npregiv(y = dat$y, z = dat$z, w = dat$w,
            nomad = TRUE, nmulti = 1L, iterate.max = 2L),
    "automatic-degree NOMAD search is not available", fixed = TRUE
  )

  dat.na <- dat
  dat.na$y[c(3L, 11L)] <- NA_real_
  formula.na <- suppressWarnings(npregiv(
    y ~ z | w, data = dat.na, na.action = na.exclude,
    nmulti = 1L, iterate.max = 2L
  ))
  expect_length(formula.na$phi, nrow(dat) - 2L)
  expect_length(fitted(formula.na), nrow(dat))
  expect_true(all(is.na(fitted(formula.na)[c(3L, 11L)])))
  expect_true(all(is.na(gradients(formula.na)[c(3L, 11L)])))
  expect_true(all(is.na(residuals(formula.na)[c(3L, 11L)])))
})

test_that("npregivderiv formula and regression-consistent accessors preserve state", {
  if (identical(iv_interface_package(), "npRmpi"))
    skip("installed MPI route sentinel owns npRmpi numerical coverage")
  dat <- iv_interface_fixture()

  native <- suppressWarnings(npregivderiv(
    y = dat$y, z = dat$z, w = dat$w,
    nmulti = 1L, iterate.max = 2L
  ))
  formula <- suppressWarnings(npregivderiv(
    y ~ z | w, data = dat,
    nmulti = 1L, iterate.max = 2L
  ))
  explicit.ll <- suppressWarnings(npregivderiv(
    y = dat$y, z = dat$z, w = dat$w, regtype = "ll",
    nmulti = 1L, iterate.max = 2L
  ))
  explicit.lc <- suppressWarnings(npregivderiv(
    y = dat$y, z = dat$z, w = dat$w, regtype = "lc",
    nmulti = 1L, iterate.max = 2L
  ))

  expect_identical(formula$phi, native$phi)
  expect_identical(formula$phi.prime, native$phi.prime)
  expect_identical(explicit.ll$phi, native$phi)
  expect_identical(explicit.ll$phi.prime, native$phi.prime)
  expect_identical(explicit.ll$bws, native$bws)
  expect_identical(formula$call[[1L]], quote(npregivderiv))
  expect_false("..." %in% names(formula$call))
  expect_null(attr(formula$call, ".Environment"))
  expect_identical(fitted(formula), formula$phi)
  expect_identical(gradients(formula), formula$phi.prime)
  expect_equal(residuals(formula), dat$y - formula$phi, tolerance = 0)
  expect_s3_class(summary(formula), "summary.npregivderiv")
  expect_identical(summary(formula)$evaluated, length(formula$norm.stop))
  expect_output(print(summary(formula)),
                "Instrumental Kernel Derivative Estimation", fixed = TRUE)
  expect_output(print(summary(formula)), "States evaluated", fixed = TRUE)
  method <- getFromNamespace("npregivderiv.default", iv_interface_package())
  expect_identical(
    eval(formals(method)$regtype),
    c("ll", "lc", "lp")
  )
  expect_null(formula$smoothing.spec$requested$regtype)
  expect_identical(formula$smoothing.spec$effective$regtype, "ll")
  expect_identical(formula$smoothing.spec$effective$degree, 1L)
  expect_identical(formula$smoothing.spec$effective$source,
                   "derivative-default")
  expect_identical(explicit.ll$smoothing.spec$requested$regtype, "ll")
  expect_identical(explicit.ll$smoothing.spec$effective$source, "explicit")
  expect_identical(explicit.lc$smoothing.spec$effective$regtype, "lc")
  expect_error(
    npregivderiv(y = dat$y, z = dat$z, w = dat$w,
                 nomad = "auto", nmulti = 1L, iterate.max = 2L),
    "NOMAD routing is not yet available", fixed = TRUE
  )
  expect_error(
    npregivderiv(y = dat$y, z = dat$z, w = dat$w, x = dat$w),
    "does not support a separate exogenous x", fixed = TRUE
  )
})

test_that("npregivderiv uses LC only for categorical-only internal stages", {
  if (identical(iv_interface_package(), "npRmpi"))
    skip("installed MPI route sentinel owns npRmpi numerical coverage")
  dat <- iv_interface_fixture(36L)
  dat$wf <- factor(ifelse(dat$w > 0, "high", "low"))

  fit <- suppressWarnings(npregivderiv(
    y = dat$y, z = dat$z, w = dat["wf"],
    nmulti = 1L, iterate.max = 2L
  ))

  expect_s3_class(fit, "npregivderiv")
  expect_identical(fit$smoothing.spec$effective$regtype, "ll")
  expect_identical(fit$stage.specs$E.y.z$regtype, "ll")
  expect_identical(fit$stage.specs$E.y.w$regtype, "lc")
  expect_identical(fit$stage.specs$residual.w$regtype, "lc")
  expect_identical(fit$stage.specs$E.y.w$degree, integer())
  expect_identical(fit$stage.specs$residual.w$degree, integer())
})
