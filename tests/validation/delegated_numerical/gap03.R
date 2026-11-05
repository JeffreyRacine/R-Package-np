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
test_that("npregiv starting values and first-state replay are coherent", {
    invisible(NULL)
    dat <- iv_forensic_fixture()
    start <- seq(-0.15, 0.15, length.out = nrow(dat))
    started <- suppressWarnings(suppressMessages(npregiv(y = dat$y, z = dat$z, w = dat$w, starting.values = start, regtype = "lc", 
        nmulti = 1L, iterate.max = 2L, stop.on.increase = FALSE)))
    expect_identical(started$starting.values.phi, start)
    expect_true(all(is.finite(started$phi.deriv.1)))
    fitted.path <- suppressWarnings(suppressMessages(npregiv(y = dat$y, z = dat$z, w = dat$w, regtype = "lc", nmulti = 1L, 
        iterate.max = 3L, stop.on.increase = FALSE, iterate.diff.tol = 0)))
    replay.bw <- fitted.path
    replay.bw$norm.index <- 1L
    grid <- data.frame(z = seq(min(dat$z), max(dat$z), length.out = 7L))
    replay <- suppressWarnings(suppressMessages(npregiv(y = dat$y, z = dat$z, w = dat$w, zeval = grid, bw = replay.bw, regtype = "lc")))
    expect_identical(replay$norm.index, 1L)
    expect_identical(replay$convergence, "BANDWIDTH_REPLAY")
    expect_identical(ncol(replay$phi.mat), 1L)
    expect_identical(length(replay$norm.stop), 1L)
    expect_identical(as.numeric(replay$phi), as.numeric(fitted.path$phi.mat[, 1L]))
    expect_length(replay$phi.eval, nrow(grid))
    grDevices::pdf()
    on.exit(grDevices::dev.off(), add = TRUE)
    expect_silent(plot(replay))
    expect_silent(plot(replay, data_overlay = TRUE, main = "IV fit", xlim = range(dat$z), lty = 2))
    expect_silent(plot(replay, gradients = TRUE))
})
