# Installed numerical contract formerly delegated from testthat.
# Run via ../run_delegated_numerical_contract.R; this is outside the default lane.
test_that("adaptive higher-order conditional objectives retain scalar oracle", {
    old_options <- options(np.messages = FALSE, np.tree = FALSE, np.largeh = FALSE, matprod = "default")
    on.exit(options(old_options), add = TRUE)
    set.seed(2026073317L)
    n <- 96L
    xdat <- data.frame(x1 = runif(n, -0.85, 0.9), x2 = runif(n, -0.8, 0.88))
    ydat <- data.frame(y = sin(1.4 * xdat$x1) + 0.3 * xdat$x2^2 + rnorm(n, sd = 0.2))
    bws <- c(31L, 29L, 27L)
    objective <- function(route, order, degree, accelerate) {
        collective(substitute(options(np.macMseries.accelerate = A),list(A=accelerate)))
        common <- list(xdat = xdat, ydat = ydat, bws = bws, bandwidth.compute = FALSE, bwtype = "adaptive_nn", regtype = "lp", 
            basis = "glp", degree = degree, bernstein.basis = FALSE, cxkertype = "gaussian", cxkerorder = order, cykertype = "gaussian", 
            cykerorder = 2L)
        if (route == "cdens_cvml") {
            state <- do.call(npcdensbw, c(common, list(bwmethod = "cv.ml")))
            return(as.numeric(collective(substitute(npRmpi:::.npcdensbw_eval_only(X, Y, S,force.local=FALSE)$objective,list(X=xdat,Y=ydat,S=state)))))
        }
        if (route == "cdens_cvls") {
            state <- do.call(npcdensbw, c(common, list(bwmethod = "cv.ls")))
            return(as.numeric(collective(substitute(npRmpi:::.npcdensbw_eval_only(X, Y, S,force.local=FALSE)$objective,list(X=xdat,Y=ydat,S=state)))))
        }
        state <- do.call(npcdistbw, c(common, list(bwmethod = "cv.ls")))
        as.numeric(collective(substitute(npRmpi:::.npcdistbw_eval_only(X,Y,bws=S,do.full.integral=TRUE,force.local=FALSE)$objective,list(X=xdat,Y=ydat,S=state))))
    }
    for (route in c("cdens_cvml", "cdens_cvls", "cdist_cvls")) {
        for (order in c(4L, 6L, 8L)) {
            for (degree in list(c(0L, 0L), c(2L, 2L))) {
                scalar <- objective(route, order, degree, FALSE)
                vector <- objective(route, order, degree, TRUE)
                expect_true(is.finite(scalar) && is.finite(vector) && max(abs(c(scalar,vector))) < .Machine$double.xmax/2)
                cat("CELL",route,order,paste(degree,collapse=","),sprintf("%.17g",scalar),sprintf("%.17g",vector),"\n", file=file.path(Sys.getenv("NP_NUMERICAL_RESULT_DIR"),"core_values.tsv"),append=TRUE)
                expect_equal(vector, scalar, tolerance = 3e-09)
            }
        }
    }
})
