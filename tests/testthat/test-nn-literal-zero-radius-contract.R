test_that("adaptive NN direct estimators reject literal zero radii locally", {
  local <- getFromNamespace(".npRmpi_with_local_regression", "npRmpi")
  x <- data.frame(x = c(0, 0, 0, 1, 2))
  y <- c(0, 1, 2, 3, 4)
  response <- data.frame(y = c(-2, -1, 0, 1, 2))

  regression_bw <- local(npregbw(
    xdat = x, ydat = y, bwtype = "adaptive_nn", bws = 2,
    bandwidth.compute = FALSE
  ))
  density_bw <- local(npudensbw(
    dat = x, bwtype = "adaptive_nn", bws = 2,
    bandwidth.compute = FALSE
  ))
  distribution_bw <- local(npudistbw(
    dat = x, bwtype = "adaptive_nn", bws = 2,
    bandwidth.compute = FALSE
  ))
  conditional_density_bw <- local(npcdensbw(
    xdat = x, ydat = response, bwtype = "adaptive_nn", bws = c(4, 2),
    bandwidth.compute = FALSE
  ))
  conditional_distribution_bw <- local(npcdistbw(
    xdat = x, ydat = response, bwtype = "adaptive_nn", bws = c(4, 2),
    bandwidth.compute = FALSE
  ))

  expect_error(local(npreg(bws = regression_bw, exdat = x)),
               "zero literal radius", fixed = TRUE)
  expect_error(local(npudens(bws = density_bw, tdat = x)),
               "zero literal radius", fixed = TRUE)
  expect_error(local(npudist(bws = distribution_bw, tdat = x)),
               "zero literal radius", fixed = TRUE)
  expect_error(
    local(npcdens(bws = conditional_density_bw, txdat = x,
                  tydat = response, exdat = x, eydat = response)),
    "zero literal radius", fixed = TRUE
  )
  expect_error(
    local(npcdist(bws = conditional_distribution_bw, txdat = x,
                  tydat = response, exdat = x, eydat = response)),
    "zero literal radius", fixed = TRUE
  )
})

test_that("beta and kernel-sum direct owners retain literal zero status", {
  local <- getFromNamespace(".npRmpi_with_local_regression", "npRmpi")
  x <- data.frame(x = c(0.1, 0.1, 0.1, 0.6, 0.9))
  y <- seq_len(nrow(x))

  for (mode in c("generalized_nn", "adaptive_nn")) {
    beta_bw <- local(npregbw(
      xdat = x, ydat = y, bwtype = mode, bws = 2,
      bandwidth.compute = FALSE, ckertype = "beta",
      ckerbound = "fixed", ckerlb = 0, ckerub = 1
    ))
    expect_error(local(npreg(bws = beta_bw, exdat = x)),
                 "zero literal radius", fixed = TRUE)
    expect_error(
      local(npksum(txdat = x, exdat = x, bws = 2, bwtype = mode,
                   ckertype = "gaussian")),
      "zero literal radius", fixed = TRUE
    )
    expect_error(
      local(npksum(txdat = x, exdat = x, bws = 2, bwtype = mode,
                   ckertype = "beta", ckerbound = "fixed",
                   ckerlb = 0, ckerub = 1)),
      "zero literal radius", fixed = TRUE
    )
  }
})

test_that("adaptive literal-zero status is rank symmetric", {
  env <- npRmpi_subprocess_env(c("NP_RMPI_NO_REUSE_SLAVES=1"))
  skip_if(is.null(env), "installed npRmpi unavailable for subprocess proof")

  for (nslaves in c(1L, 3L)) {
    marker <- sprintf("NPRMPI_ADAPTIVE_ZERO_RADIUS_OK_%d", nslaves)
    result <- npRmpi_run_rscript_subprocess(
      lines = c(
        "suppressPackageStartupMessages(library(npRmpi))",
        "options(npRmpi.autodispatch=TRUE, np.messages=FALSE, np.tree=FALSE)",
        sprintf("npRmpi.init(nslaves=%dL, quiet=TRUE)", nslaves),
        "on.exit(try(npRmpi.quit(force=TRUE), silent=TRUE), add=TRUE)",
        "x <- data.frame(x=c(0,0,0,1,2)); y <- 0:4",
        "bad <- npregbw(xdat=x,ydat=y,bwtype='adaptive_nn',bws=2,bandwidth.compute=FALSE)",
        "good <- npregbw(xdat=x,ydat=y,bwtype='adaptive_nn',bws=3,bandwidth.compute=FALSE)",
        "condition <- tryCatch(npreg(bws=bad,exdat=x), error=identity)",
        "stopifnot(inherits(condition,'np_nn_zero_radius'), identical(condition$variable,'x'), identical(condition$lookup.k,2L), identical(condition$excluded,1L), grepl('zero literal radius',conditionMessage(condition),fixed=TRUE))",
        "fit <- npreg(bws=good,exdat=x)",
        "stopifnot(length(fit$mean)==nrow(x),all(is.finite(fit$mean)))",
        "smoke <- npRmpi:::.npRmpi_spmd_tiny_smoke('after-adaptive-radius-error',comm=1L)",
        "stopifnot(is.list(smoke),isTRUE(smoke$ok))",
        "z <- data.frame(group=factor(rep(c('a','b'),6L)),good=seq_len(12L)/12,tied=c(rep(0,4L),1:8)); response <- 1+sin(seq_len(12L))/4+z$good",
        "gnn <- npregbw(xdat=z,ydat=response,bwtype='generalized_nn',bws=c(.2,8,4),regtype='lc',bandwidth.compute=FALSE)",
        "gnn.fit <- npreg(bws=gnn,txdat=z,tydat=response,se=FALSE)",
        "set.seed(527); seed <- .Random.seed; sig.condition <- tryCatch(npsigtest(gnn.fit,B=9L,joint=FALSE),error=identity)",
        "stopifnot(inherits(sig.condition,'np_nn_zero_radius'),identical(sig.condition$variable,'tied'),grepl('bootstrap geometry preflight',sig.condition$stage,fixed=TRUE),identical(.Random.seed,seed))",
        "smoke <- npRmpi:::.npRmpi_spmd_tiny_smoke('after-npsig-radius-error',comm=1L)",
        "stopifnot(is.list(smoke),isTRUE(smoke$ok))",
        sprintf("cat('%s\\n')", marker)
      ),
      timeout = 90L,
      env = env
    )
    info <- paste(result$output, collapse = "\n")
    expect_equal(result$status, 0L, info = info)
    expect_true(any(grepl(marker, result$output, fixed = TRUE)), info = info)
  }
})
