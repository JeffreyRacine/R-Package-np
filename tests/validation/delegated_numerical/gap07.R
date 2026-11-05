# Installed numerical contract formerly delegated from testthat.
# Run via ../run_delegated_numerical_contract.R; this is outside the default lane.
iv_interface_fixture <- function(n = 24L) {
    set.seed(9010)
    w <- rnorm(n)
    v <- rnorm(n, sd = 0.24)
    z <- 0.4 * w + v
    y <- z^2 - 0.45 * v + rnorm(n, sd = 0.06)
    data.frame(y = y, z = z, w = w)
}
iv_interface_package <- function() {
    if ("npRmpi" %in% loadedNamespaces()) 
        "npRmpi"
    else "npRmpi"
}
test_that("npregiv formula, modern degree spelling, and accessors preserve state", {
    invisible(NULL)
    dat <- iv_interface_fixture()
    native <- suppressWarnings(npregiv(y = dat$y, z = dat$z, w = dat$w, nmulti = 1L, iterate.max = 2L))
    formula <- suppressWarnings(npregiv(y ~ z | w, data = dat, nmulti = 1L, iterate.max = 2L))
    modern.ll <- suppressWarnings(npregiv(y = dat$y, z = dat$z, w = dat$w, regtype = "ll", nmulti = 1L, iterate.max = 2L))
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
    expect_error(npregiv(y = dat$y, z = dat$z, w = dat$w, p = 1L, regtype = "lc", nmulti = 1L, iterate.max = 2L), "conflicting smoothing controls", 
        fixed = TRUE)
    expect_error(npregiv(y = dat$y, z = dat$z, w = dat$w, nomad = TRUE, nmulti = 1L, iterate.max = 2L), "automatic-degree NOMAD search is not available", 
        fixed = TRUE)
    dat.na <- dat
    dat.na$y[c(3L, 11L)] <- NA_real_
    formula.na <- suppressWarnings(npregiv(y ~ z | w, data = dat.na, na.action = na.exclude, nmulti = 1L, iterate.max = 2L))
    expect_length(formula.na$phi, nrow(dat) - 2L)
    expect_length(fitted(formula.na), nrow(dat))
    expect_true(all(is.na(fitted(formula.na)[c(3L, 11L)])))
    expect_true(all(is.na(gradients(formula.na)[c(3L, 11L)])))
    expect_true(all(is.na(residuals(formula.na)[c(3L, 11L)])))
})
