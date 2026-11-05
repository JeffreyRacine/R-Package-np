# Installed numerical contract formerly delegated from testthat.
# Run via ../run_delegated_numerical_contract.R; this is outside the default lane.
iv_forensic_fixture <- function(n = 26L, seed = 20260721L) {
    set.seed(seed)
    w <- rnorm(n)
    v <- rnorm(n, sd = 0.22)
    z <- 0.45 * w + v
    y <- z^2 - 0.4 * v + rnorm(n, sd = 0.055)
    data.frame(y, z, w)
}
iv_forensic_skip_mpi_numerics <- function() {
    if ("npRmpi" %in% loadedNamespaces()) 
        skip("installed one- and multi-worker sentinels own npRmpi numerical coverage")
}
test_that("npregivderiv evaluation is output-only and S3 accessors stay training-aligned", {
    invisible(NULL)
    dat <- iv_forensic_fixture()
    args <- list(nmulti = 1L, iterate.max = 3L, iterate.break = FALSE, stop.on.increase = FALSE, random.seed = 20260721L)
    set.seed(73)
    caller.seed <- .Random.seed
    training <- suppressWarnings(suppressMessages(do.call(npregivderiv, c(list(y = dat$y, z = dat$z, w = dat$w), args))))
    expect_identical(.Random.seed, caller.seed)
    reversed <- suppressWarnings(suppressMessages(do.call(npregivderiv, c(list(y = dat$y, z = dat$z, w = dat$w, zeval = rev(dat$z)), 
        args))))
    grid <- data.frame(z = seq(min(dat$z), max(dat$z), length.out = 9L))
    evaluated <- suppressWarnings(suppressMessages(do.call(npregivderiv, c(list(y = dat$y, z = dat$z, w = dat$w, zeval = grid), 
        args))))
    formula.evaluated <- suppressWarnings(suppressMessages(do.call(npregivderiv, c(list(y ~ z | w, data = dat, newdata = grid), 
        args))))
    for (field in c("phi", "phi.prime", "phi.mat", "phi.prime.mat", "norm.stop", "num.iterations")) {
        expect_identical(reversed[[field]], training[[field]])
        expect_identical(evaluated[[field]], training[[field]])
    }
    expect_identical(reversed$phi.prime.eval, rev(training$phi.prime))
    expect_identical(formula.evaluated$phi.prime.eval, evaluated$phi.prime.eval)
    expect_identical(fitted(evaluated), evaluated$phi)
    expect_identical(gradients(evaluated), evaluated$phi.prime)
    expect_equal(residuals(evaluated), dat$y - evaluated$phi, tolerance = 0)
    expect_error(do.call(npregivderiv, c(list(y = dat$y, z = dat$z, w = dat$w, weval = rev(dat$w)), args)), "weval cannot differ", 
        fixed = TRUE)
    grDevices::pdf()
    on.exit(grDevices::dev.off(), add = TRUE)
    expect_silent(plot(evaluated))
    expect_silent(plot(evaluated, gradients = FALSE))
    expect_silent(plot(evaluated, data_overlay = TRUE, gradients = FALSE, main = "IV derivative fit", lty = 2))
    expect_output(print(summary(evaluated)), "Bandwidth for E(y-phi(z)|w):", fixed = TRUE)
})
