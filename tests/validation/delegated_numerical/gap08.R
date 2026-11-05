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
test_that("npregivderiv formula and regression-consistent accessors preserve state", {
    invisible(NULL)
    dat <- iv_interface_fixture()
    native <- suppressWarnings(npregivderiv(y = dat$y, z = dat$z, w = dat$w, nmulti = 1L, iterate.max = 2L))
    formula <- suppressWarnings(npregivderiv(y ~ z | w, data = dat, nmulti = 1L, iterate.max = 2L))
    explicit.ll <- suppressWarnings(npregivderiv(y = dat$y, z = dat$z, w = dat$w, regtype = "ll", nmulti = 1L, iterate.max = 2L))
    explicit.lc <- suppressWarnings(npregivderiv(y = dat$y, z = dat$z, w = dat$w, regtype = "lc", nmulti = 1L, iterate.max = 2L))
    expect_identical(formula$phi, native$phi)
    expect_identical(formula$phi.prime, native$phi.prime)
    expect_identical(explicit.ll$phi, native$phi)
    expect_identical(explicit.ll$phi.prime, native$phi.prime)
    volatile <- c("wall_elapsed_sec", "compute_elapsed_sec", "user_sec", "system_sec", "timestamp_start", "timestamp_end")
    clean <- function(bws) {
        profile <- bws[["timing.profile", exact = TRUE]]
        if (!is.null(profile)) {
            for (name in intersect(volatile, names(profile))) {
                expect_type(profile[[name]], "double")
                expect_length(profile[[name]], 1L)
                expect_true(is.finite(profile[[name]]))
                profile[[name]] <- 0
            }
            bws$timing.profile <- profile
        }
        bws
    }
    expect_identical(clean(explicit.ll$bws), clean(native$bws))
    expect_identical(formula$call[[1L]], quote(npregivderiv))
    expect_false("..." %in% names(formula$call))
    expect_null(attr(formula$call, ".Environment"))
    expect_identical(fitted(formula), formula$phi)
    expect_identical(gradients(formula), formula$phi.prime)
    expect_equal(residuals(formula), dat$y - formula$phi, tolerance = 0)
    expect_s3_class(summary(formula), "summary.npregivderiv")
    expect_identical(summary(formula)$evaluated, length(formula$norm.stop))
    expect_output(print(summary(formula)), "Instrumental Kernel Derivative Estimation", fixed = TRUE)
    expect_output(print(summary(formula)), "States evaluated", fixed = TRUE)
    method <- getFromNamespace("npregivderiv.default", iv_interface_package())
    expect_identical(eval(formals(method)$regtype), c("ll", "lc", "lp"))
    expect_null(formula$smoothing.spec$requested$regtype)
    expect_identical(formula$smoothing.spec$effective$regtype, "ll")
    expect_identical(formula$smoothing.spec$effective$degree, 1L)
    expect_identical(formula$smoothing.spec$effective$source, "derivative-default")
    expect_identical(explicit.ll$smoothing.spec$requested$regtype, "ll")
    expect_identical(explicit.ll$smoothing.spec$effective$source, "explicit")
    expect_identical(explicit.lc$smoothing.spec$effective$regtype, "lc")
    expect_error(npregivderiv(y = dat$y, z = dat$z, w = dat$w, nomad = "auto", nmulti = 1L, iterate.max = 2L), "NOMAD routing is not yet available", 
        fixed = TRUE)
    expect_error(npregivderiv(y = dat$y, z = dat$z, w = dat$w, x = dat$w), "does not support a separate exogenous x", fixed = TRUE)
})
