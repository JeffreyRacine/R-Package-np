test_that("fixed tree ownership cannot depend on sample size", {
  header <- paste(
    readLines(npRmpi_test_source_path("src", "tree_capability.h"),
              warn = FALSE),
    collapse = "\n"
  )
  source <- paste(
    readLines(npRmpi_test_source_path("src", "jksum.c"), warn = FALSE),
    collapse = "\n"
  )
  selector <- sub(
    "(?s).*?(fixed_tree_common_geometry[[:space:]]*=[[:space:]]*np_ks_tree_use.*?const char \\* const owner_oracle).*",
    "\\1", source, perl = TRUE
  )

  expect_false(grepl(
    "\\b(num_obs|nobs|sample_size|sample_count)\\b",
    header, perl = TRUE
  ))
  expect_false(grepl(
    "num_obs_(train|eval)[[:space:]]*(<=|>=|<|>|==|!=)",
    selector, perl = TRUE
  ))
})

test_that("exact inactive MPI tree dimensions do not trigger K(0) arithmetic", {
  source <- paste(
    readLines(npRmpi_test_source_path("src", "jksum.c"), warn = FALSE),
    collapse = "\n"
  )
  active_plan <- sub(
    "(?s).*?(const int tree_use_active_dims[[:space:]]*=.*?for \\(ii = 0; ii < p_nvar; ii\\+\\+\\)).*",
    "\\1", source, perl = TRUE
  )
  kernel_plan <- sub(
    "(?s).*?(const int use_largeh = any_cont_largeh.*?np_ckernelv\\().*",
    "\\1", source, perl = TRUE
  )

  expect_match(active_plan, "fixed_tree_nonpruning_dimension",
               fixed = TRUE)
  expect_match(active_plan, "cont_largeh_active", fixed = TRUE)
  expect_false(grepl("fixed_tree_nonpruning_dimension", kernel_plan,
                     fixed = TRUE))
})

test_that("MPI fixed tree capability is rank-symmetric", {
  skip_on_cran()
  env <- npRmpi_subprocess_env(c("NP_RMPI_NO_REUSE_SLAVES=1"))
  skip_if(is.null(env), "local npRmpi install unavailable for subprocess smoke")

  for (nslaves in c(1L, 3L)) {
    marker <- paste0("TREE_CAPABILITY_MPI_OWNER_OK_", nslaves)
    result <- npRmpi_run_rscript_subprocess(
      lines = c(
        "suppressPackageStartupMessages(library(npRmpi))",
        sprintf("npRmpi.init(nslaves=%dL, quiet=TRUE)", nslaves),
        "options(np.messages=FALSE, np.largeh=FALSE, np.macMseries.accelerate=FALSE)",
        "set.seed(20260820L)",
        "n <- 61L",
        "x <- data.frame(x1=runif(n),x2=runif(n))",
        "y <- data.frame(y=0.4+x$x1-0.6*x$x2+rnorm(n,sd=0.2))",
        "responses <- cbind(1,y$y)",
        "weights <- cbind(1,seq_len(n)/n)",
        "xr <- vapply(x,function(value) diff(range(value)),numeric(1L))",
        "yr <- diff(range(y$y))",
        "full <- 1.25*xr/sqrt(5)",
        "prunable <- 0.20*xr/sqrt(5)",
        "mixed <- c(full[[1L]],prunable[[2L]])",
        "run_kws <- function(tree,case) {",
        "  options(np.tree=tree)",
        "  do.call(npksum,c(list(txdat=x,tydat=responses,weights=weights,ckertype='epanechnikov',leave.one.out=TRUE,return.kernel.weights=TRUE),case))",
        "}",
        "cases <- list(",
        "  normal_full=list(bws=full,operator='normal'),",
        "  normal_prunable=list(bws=prunable,operator='normal'),",
        "  normal_mixed=list(bws=mixed,operator='normal'),",
        "  convolution_full=list(bws=full,operator='convolution'),",
        "  derivative_full=list(bws=full,operator='derivative'),",
        "  integral_full=list(bws=full,operator='integral'),",
        "  union_prunable=list(bws=0.30*xr,operator='convolution',permutation.operator='derivative',return.derivative.kernel.weights=TRUE)",
        ")",
        "for (case in cases) {",
        "  dense <- run_kws(FALSE,case); tree <- run_kws(TRUE,case)",
        "  for (field in c('ksum','kw','p.ksum','p.kw'))",
        "    stopifnot(isTRUE(all.equal(tree[[field]],dense[[field]],tolerance=2e-11)))",
        "}",
        "owner <- function(case) {",
        "  Sys.setenv(NP_KWS_TREE_OWNER_ORACLE='1',NP_KWS_TREE_DIMENSION_ORACLE='1')",
        "  messages <- capture.output(run_kws(TRUE,case),type='message')",
        "  Sys.unsetenv(c('NP_KWS_TREE_OWNER_ORACLE','NP_KWS_TREE_DIMENSION_ORACLE'))",
        "  messages",
        "}",
        "stopifnot(any(grepl('common=1 capability=1 owner=dense',owner(cases$normal_full),fixed=TRUE)))",
        "stopifnot(any(grepl('common=1 capability=2 owner=tree',owner(cases$normal_prunable),fixed=TRUE)))",
        "stopifnot(any(grepl('dimensions=2 inactive=1 active=1',owner(cases$normal_mixed),fixed=TRUE)))",
        "stopifnot(any(grepl('common=1 capability=2 owner=tree',owner(cases$union_prunable),fixed=TRUE)))",
        "core <- getFromNamespace('.npregbw_call_fixed_degree_core','npRmpi')",
        "reg_owner <- function(degree,bandwidth) {",
        "  bw <- npregbw(xdat=x,ydat=y$y,bws=bandwidth,bandwidth.compute=FALSE,bwtype='fixed',bwmethod='cv.ls',ckertype='epanechnikov',regtype='lp',degree=degree)",
        "  options(np.tree=TRUE); Sys.setenv(NP_LP_TREE_ORACLE='1')",
        "  messages <- capture.output(tree <- core(xdat=x,ydat=y$y,bws=bw,eval.only=TRUE),type='message')",
        "  Sys.unsetenv('NP_LP_TREE_ORACLE'); options(np.tree=FALSE)",
        "  dense <- core(xdat=x,ydat=y$y,bws=bw,eval.only=TRUE)",
        "  stopifnot(isTRUE(all.equal(tree$objective,dense$objective,tolerance=2e-11)))",
        "  any(grepl('NP_LP_TREE_ORACLE',messages,fixed=TRUE))",
        "}",
        "stopifnot(!reg_owner(c(1L,1L),full))",
        "stopifnot(reg_owner(c(1L,1L),prunable))",
        "stopifnot(reg_owner(c(1L,0L),full))",
        "reg_bw <- npregbw(xdat=x,ydat=y$y,bws=mixed,bandwidth.compute=FALSE,bwtype='fixed',bwmethod='cv.ls',ckertype='epanechnikov',regtype='lc')",
        "udens_bw <- npudensbw(dat=x,bws=mixed,bandwidth.compute=FALSE,bwmethod='cv.ml',ckertype='epanechnikov')",
        "udist_bw <- npudistbw(dat=x,bws=mixed,bandwidth.compute=FALSE,bwmethod='cv.cdf',ckertype='epanechnikov')",
        "cdens_bw <- npcdensbw(xdat=x,ydat=y,bws=c(1.25*yr/sqrt(5),full),bandwidth.compute=FALSE,bwmethod='cv.ml',regtype='lc',cxkertype='epanechnikov',cykertype='epanechnikov')",
        "public <- function(tree) {",
        "  options(np.tree=tree)",
        "  list(reg=npreg(bws=reg_bw,se=TRUE,gradients=FALSE),udens_objective=npudensbw(dat=x,bws=udens_bw,bandwidth.compute=TRUE,eval.only=TRUE)$fval,udens=npudens(bws=udens_bw),udist=npudist(bws=udist_bw),cdens=npRmpi:::.npcdensbw_eval_only(x,y,cdens_bw)$objective,fit=npcdens(bws=cdens_bw,txdat=x,tydat=y,se=TRUE,gradients=TRUE))",
        "}",
        "dense_public <- public(FALSE); tree_public <- public(TRUE)",
        "for (field in c('mean','merr','resid')) stopifnot(isTRUE(all.equal(tree_public$reg[[field]],dense_public$reg[[field]],tolerance=2e-11)))",
        "stopifnot(isTRUE(all.equal(tree_public$udens_objective,dense_public$udens_objective,tolerance=2e-11)))",
        "for (field in c('dens','derr','log_likelihood')) stopifnot(isTRUE(all.equal(tree_public$udens[[field]],dense_public$udens[[field]],tolerance=2e-11)))",
        "for (field in c('dist','derr')) stopifnot(isTRUE(all.equal(tree_public$udist[[field]],dense_public$udist[[field]],tolerance=2e-11)))",
        "stopifnot(isTRUE(all.equal(tree_public$cdens,dense_public$cdens,tolerance=2e-11)))",
        "for (field in c('condens','conderr','congrad','congerr'))",
        "  stopifnot(isTRUE(all.equal(tree_public$fit[[field]],dense_public$fit[[field]],tolerance=2e-11)))",
        "npRmpi.quit(force=TRUE)",
        sprintf("cat('%s\\n')", marker)
      ),
      timeout = 180L,
      env = env
    )
    expect_equal(result$status, 0L,
                 info = paste(result$output, collapse = "\n"))
    expect_true(any(grepl(marker, result$output, fixed = TRUE)),
                info = paste(result$output, collapse = "\n"))
  }
})
