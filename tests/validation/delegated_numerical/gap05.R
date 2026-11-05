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
test_that("npregiv owns categorical kernels and bandwidth scaling arguments", {
    invisible(NULL)
    dat <- iv_forensic_fixture(34L)
    wf <- factor(ifelse(dat$w > 0, "high", "low"))
    fit <- suppressWarnings(suppressMessages(npregiv(y = dat$y, z = dat$z, w = data.frame(wf), regtype = "lc", nmulti = 1L, 
        iterate.max = 2L, ukertype = "aitchisonaitken", bandwidth.divide = FALSE)))
    upper <- getFromNamespace("uMaxL", "npRmpi")(nlevels(wf), "aitchisonaitken")
    expect_s3_class(fit, "npregiv")
    expect_lte(as.numeric(fit$bw.E.y.w[[1L]]), upper)
    expect_identical(fit$stage.specs$w$regtype, "lc")
})
