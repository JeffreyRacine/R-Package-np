# Installed numerical contract formerly delegated from testthat.
# Run via ../run_delegated_numerical_contract.R; this is outside the default lane.
test_that("accelerated Gaussian CVCDF remains numerically tied to scalar CVCDF", {
    owner <- environmentName(environment(npudistbw))
    invisible(NULL)
    set.seed(20260727L)
    dat <- data.frame(x1 = rnorm(300L), x2 = rnorm(300L))
    bw <- npudistbw(dat = dat, bws = c(0.35, 0.38), bandwidth.compute = FALSE, bwmethod = "cv.cdf")
    evaluate <- function(accelerate) {
      collective(substitute({
        options(np.macMseries.accelerate=A)
        getFromNamespace("npudistbw.dbandwidth","npRmpi")(dat=D,bws=S,
          bandwidth.compute=TRUE,eval.only=TRUE,nmulti=1L,powell.remin=FALSE)$fval[[1L]]
      }, list(A=accelerate,D=dat,S=bw)))
    }
    scalar <- evaluate(FALSE)
    accelerated <- evaluate(TRUE)
    cat("CELL cvcdf",sprintf("%.17g",scalar),sprintf("%.17g",accelerated),"\n",file=file.path(Sys.getenv("NP_NUMERICAL_RESULT_DIR"),"core_values.tsv"))
    expect_true(is.finite(scalar))
    expect_true(is.finite(accelerated))
    expect_equal(accelerated, scalar, tolerance = 1e-10)
})
