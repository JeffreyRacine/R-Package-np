test_that("fixed full-support tree candidates preserve canonical MPI LP objectives", {
  skip_on_cran()
  env <- npRmpi_subprocess_env(c("NP_RMPI_NO_REUSE_SLAVES=1"))
  skip_if(is.null(env), "local npRmpi install unavailable for subprocess smoke")

  marker <- "REGRESSION_TREE_FULL_SUPPORT_OWNER_OK"
  res <- npRmpi_run_rscript_subprocess(
    lines = c(
      "suppressPackageStartupMessages(library(npRmpi))",
      "npRmpi.init(nslaves = 1L, quiet = TRUE)",
      "options(np.messages = FALSE, np.largeh = TRUE, np.macMseries.accelerate = FALSE)",
      "set.seed(20260812L)",
      "n <- 89L",
      "all_x <- data.frame(x1=runif(n,-1,1), x2=runif(n,-1,1), x3=runif(n,-1,1))",
      "y <- 0.5 + 0.7*all_x$x1 - 0.4*all_x$x2 + 0.3*all_x$x1*all_x$x2 + rnorm(n,sd=0.25)",
      "core <- getFromNamespace('.npregbw_call_fixed_degree_core','npRmpi')",
      "check_case <- function(p, degree, bernstein, kernel, order) {",
      "  x <- all_x[seq_len(p)]",
      "  range_x <- vapply(x, function(value) diff(range(value)), numeric(1L))",
      "  bw_args <- list(xdat=x, ydat=y, regtype='lp', degree=degree, bernstein.basis=bernstein, bwmethod='cv.ls', bwtype='fixed', ckertype=kernel, bandwidth.compute=FALSE, bws=10*range_x)",
      "  if (kernel != 'uniform') bw_args$ckerorder <- order",
      "  bw <- do.call(npregbw, bw_args)",
      "  options(np.tree=FALSE)",
      "  dense <- core(xdat=x, ydat=y, bws=bw, eval.only=TRUE)",
      "  options(np.tree=TRUE)",
      "  tree <- core(xdat=x, ydat=y, bws=bw, eval.only=TRUE)",
      "  stopifnot(isTRUE(all.equal(tree$objective,dense$objective,tolerance=2e-11)))",
      "  stopifnot(identical(tree$num.feval.fast,dense$num.feval.fast))",
      "}",
      "check_case(1L,0L,FALSE,'epanechnikov',2L)",
      "check_case(2L,c(1L,1L),TRUE,'uniform',2L)",
      "check_case(2L,c(2L,2L),FALSE,'epanechnikov',4L)",
      "check_case(3L,rep(2L,3L),TRUE,'epanechnikov',6L)",
      "check_case(3L,rep(3L,3L),TRUE,'epanechnikov',8L)",
      "mixed_x <- data.frame(x1=all_x$x1, unordered=factor(rep(c('a','b','c'),length.out=n)), ordered=ordered(rep(1:4,length.out=n)))",
      "mixed_bw <- npregbw(xdat=mixed_x, ydat=y, regtype='lp', degree=1L, bernstein.basis=TRUE, bwmethod='cv.ls', bwtype='fixed', ckertype='epanechnikov', bandwidth.compute=FALSE, bws=c(10*diff(range(mixed_x$x1)),0.3,0.3))",
      "options(np.tree=FALSE)",
      "mixed_dense <- core(xdat=mixed_x, ydat=y, bws=mixed_bw, eval.only=TRUE)",
      "options(np.tree=TRUE)",
      "mixed_tree <- core(xdat=mixed_x, ydat=y, bws=mixed_bw, eval.only=TRUE)",
      "stopifnot(isTRUE(all.equal(mixed_tree$objective,mixed_dense$objective,tolerance=2e-11)))",
      "stopifnot(identical(mixed_tree$num.feval.fast,mixed_dense$num.feval.fast))",
      "method_x <- all_x[1:2]",
      "method_y <- 0.4 + 0.7*method_x$x1 - 0.5*method_x$x2 + rnorm(n,sd=0.3)",
      "method_h <- 10*vapply(method_x,function(value) diff(range(value)),numeric(1L))",
      "check_method <- function(method, objective='ls') {",
      "  bw <- npregbw(xdat=method_x, ydat=method_y, regtype='lp', degree=c(1L,1L), bernstein.basis=TRUE, bwmethod=method, bwtype='fixed', ckertype='epanechnikov', bandwidth.compute=FALSE, bws=method_h)",
      "  options(np.tree=FALSE)",
      "  dense <- core(xdat=method_x, ydat=method_y, bws=bw, eval.only=TRUE, objective=objective)",
      "  options(np.tree=TRUE)",
      "  tree <- core(xdat=method_x, ydat=method_y, bws=bw, eval.only=TRUE, objective=objective)",
      "  stopifnot(isTRUE(all.equal(tree$objective,dense$objective,tolerance=2e-11)))",
      "  stopifnot(identical(tree$num.feval.fast,dense$num.feval.fast))",
      "}",
      "check_method('cv.ls','ls')",
      "check_method('cv.ls','ks')",
      "check_method('cv.aic','ls')",
      "lsq_h <- 10*vapply(method_x,function(value) diff(range(value)),numeric(1L))",
      "lsq_args <- list(xdat=method_x, ydat=method_y, scale=rep.int(1,n), tau=0.5, bws=lsq_h, regtype='lp', degree=c(2L,2L), bernstein.basis=TRUE, ckertype='epanechnikov', bandwidth.compute=FALSE)",
      "options(np.tree=FALSE)",
      "lsq_dense <- do.call(nplsqregbw,lsq_args)",
      "options(np.tree=TRUE)",
      "lsq_tree <- do.call(nplsqregbw,lsq_args)",
      "stopifnot(isTRUE(all.equal(lsq_tree$objective,lsq_dense$objective,tolerance=2e-11)))",
      "stopifnot(identical(lsq_tree$num.feval.fast,lsq_dense$num.feval.fast))",
      sprintf("cat('%s\\n')", marker)
    ),
    timeout = 180L,
    env = env
  )

  expect_equal(res$status, 0L, info = paste(res$output, collapse = "\n"))
  expect_true(any(grepl(marker, res$output, fixed = TRUE)),
              info = paste(res$output, collapse = "\n"))
})
