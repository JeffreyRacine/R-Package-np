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
test_that("multivariate npregiv derivatives are coordinate-specific", {
    invisible(NULL)
    set.seed(812)
    n <- 38L
    w1 <- rnorm(n)
    w2 <- rnorm(n)
    v1 <- rnorm(n, sd = 0.12)
    v2 <- rnorm(n, sd = 0.12)
    z1 <- 0.7 * w1 + 0.1 * w2 + v1
    z2 <- -0.1 * w1 + 0.7 * w2 + v2
    y <- z1^2 + 0.5 * z2^2 - 0.25 * v1 + 0.1 * v2 + rnorm(n, sd = 0.04)
    z <- data.frame(z1, z2)
    w <- data.frame(w1, w2)
    orders <- getFromNamespace(".np_iv_partial_orders", "npRmpi")
    expect_identical(orders(2L, 1L), diag(1L, 2L, 2L))
    expect_identical(orders(2L, 2L), diag(2L, 2L, 2L))
    ll <- suppressWarnings(suppressMessages(npregiv(y = y, z = z, w = w, zeval = data.frame(z1 = 0, z2 = 0), regtype = "ll", 
        nmulti = 1L, iterate.max = 2L, stop.on.increase = FALSE)))
    expect_identical(dim(ll$phi.deriv.1), c(n, 2L))
    expect_identical(colnames(ll$phi.deriv.1), names(z))
    expect_identical(dim(ll$phi.deriv.eval.1), c(1L, 2L))
    expect_identical(colnames(ll$phi.deriv.eval.1), names(z))
    expect_true(all(is.finite(ll$phi.deriv.1)))
    tikh <- suppressWarnings(suppressMessages(npregiv(y = y, z = z, w = w, method = "Tikhonov", regtype = "lc", alpha = 0.02, 
        alpha.iter = 0.02, nmulti = 1L, return.weights.phi = TRUE, return.weights.phi.deriv.1 = TRUE)))
    expect_identical(dim(tikh$phi.deriv.1), c(n, 2L))
    expect_identical(colnames(tikh$phi.deriv.1), names(z))
    expect_named(tikh$phi.deriv.1.weights, names(z))
    for (j in seq_along(z)) {
        expect_equal(as.numeric(tikh$phi.deriv.1.weights[[j]] %*% y), tikh$phi.deriv.1[, j], tolerance = 1e-11)
    }
})
