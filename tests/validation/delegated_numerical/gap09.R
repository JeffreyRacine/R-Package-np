# Installed numerical contract formerly delegated from testthat.
# Run via ../run_delegated_numerical_contract.R; this is outside the default lane.
iv_summary_namespace <- function() {
    if ("npRmpi" %in% loadedNamespaces()) 
        "npRmpi"
    else "npRmpi"
}
iv_summary_formatter <- function() {
    getFromNamespace(".np_iv_regression_summary_text", iv_summary_namespace())
}
test_that("public IV summaries replace only the smoothing block", {
    invisible(NULL)
    set.seed(2718)
    n <- 22L
    w <- rnorm(n)
    v <- rnorm(n, sd = 0.18)
    z <- 0.55 * w + v
    y <- sin(z) - 0.3 * v + rnorm(n, sd = 0.04)
    args <- list(y = y, z = z, w = w, nmulti = 1L, iterate.max = 2L, stop.on.increase = FALSE, random.seed = 2718L)
    ordinary <- suppressWarnings(suppressMessages(do.call(npregiv, c(args, list(regtype = "ll")))))
    derivative <- suppressWarnings(suppressMessages(do.call(npregivderiv, c(args, list(iterate.break = FALSE, regtype = "ll")))))
    for (fit in list(ordinary, derivative)) {
        object <- summary(fit)
        output <- capture.output(visible <- withVisible(print(object)))
        expect_false(visible$visible)
        expect_identical(visible$value, object)
        expect_identical(grep("^Kernel Regression Estimator:", output, value = TRUE), "Kernel Regression Estimator: Local-Linear")
        expect_false(any(grepl("^Local polynomial order \\(p\\):", output)))
        expect_false(any(grepl("^Regression smoothing type:", output)))
        expect_false(any(grepl("^Regression smoothing degree:", output)))
        expect_true(any(grepl("^Training observations:", output)))
        expect_true(any(grepl("^Regularization method:", output)))
    }
    malformed <- summary(derivative)
    malformed$smoothing.spec <- list(effective = list(regtype = "ll", degree = 2L))
    malformed.output <- capture.output(print(malformed))
    expect_false(any(grepl("^Kernel Regression Estimator:", malformed.output)))
    expect_true(any(grepl("^Training observations:", malformed.output)))
    expect_true(any(grepl("^Regularization method:", malformed.output)))
})
