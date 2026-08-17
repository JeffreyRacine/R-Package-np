test_that("conditional generalized-NN CVML uses exact mapped-occurrence geometry", {
  env <- npRmpi_subprocess_env(c("NP_RMPI_NO_REUSE_SLAVES=1"))
  skip_if(is.null(env), "installed npRmpi unavailable for subprocess proof")
  nslaves <- as.integer(Sys.getenv("NP_RMPI_GNN_GEOMETRY_TEST_NSLAVES", "1"))
  stopifnot(nslaves %in% c(1L, 3L))
  ok_tag <- "NPRMPI_CONDITIONAL_GNN_CVML_GEOMETRY_OK"
  lines <- c(
    "suppressPackageStartupMessages(library(npRmpi))",
    "options(npRmpi.autodispatch=TRUE, npRmpi.reuse.slaves=FALSE, np.messages=FALSE, np.tree=FALSE)",
    sprintf("npRmpi.init(nslaves=%dL, quiet=TRUE)", nslaves),
    "on.exit(try(npRmpi.quit(force=TRUE), silent=TRUE), add=TRUE)",
    "radius <- function(train,evaluation,k,excluded) sort(abs(train[setdiff(seq_along(train),excluded)]-evaluation), method='radix')[[k]]",
    "oracle <- function(x,y,kx,ky,degree=0L) {",
    "  x <- as.matrix(x); y <- as.numeric(y); n <- nrow(x)",
    "  fitted <- vapply(seq_len(n), function(i) {",
    "    hx <- vapply(seq_len(ncol(x)), function(l) radius(x[,l],x[i,l],kx[[l]],i), numeric(1L))",
    "    hy <- radius(y,y[[i]],ky,i)",
    "    wx <- rep.int(1,n)",
    "    for (l in seq_len(ncol(x))) wx <- wx*dnorm((x[i,l]-x[,l])/hx[[l]])/hx[[l]]",
    "    wy <- dnorm((y[[i]]-y)/hy)/hy; wx[[i]] <- 0; wy[[i]] <- 0",
    "    if (degree==0L) influence <- wx/sum(wx) else {",
    "      design <- cbind(1,x[,1L]-x[i,1L]); gram <- crossprod(design,wx*design)",
    "      influence <- drop(c(1,0)%*%solve(gram,t(design)*rep(wx,each=2L)))",
    "    }",
    "    sum(influence*wy)",
    "  }, numeric(1L))",
    "  stopifnot(all(is.finite(fitted)),all(fitted>0)); sum(log(fitted))",
    "}",
    "make_bw <- function(x,y,kx,ky,degree) {",
    "  args <- list(xdat=as.data.frame(x),ydat=data.frame(y=y),bws=c(ky,kx),bandwidth.compute=FALSE,bwmethod='cv.ml',bwtype='generalized_nn',bwscaling=FALSE,cxkertype='gaussian',cykertype='gaussian',regtype=if(degree==0L)'lc' else 'lp')",
    "  if (degree>0L) {args$basis <- 'glp'; args$degree <- rep.int(degree,ncol(as.matrix(x)))}",
    "  do.call(npcdensbw,args)",
    "}",
    "set.seed(20260816); n <- 13L",
    "x1 <- sort(runif(n,-0.9,1.1)); x2 <- runif(n,-0.6,0.8)",
    "y <- 0.8*x1-0.35*x2+rnorm(n,sd=0.17)",
    "cases <- list(lc1=list(x=data.frame(x1=x1),kx=4L,ky=5L,degree=0L),lc2=list(x=data.frame(x1=x1,x2=x2),kx=c(4L,6L),ky=5L,degree=0L),ll1=list(x=data.frame(x1=x1),kx=6L,ky=5L,degree=1L))",
    "for (case_name in names(cases)) {",
    "  case <- cases[[case_name]]",
    "  expected <- oracle(case$x,y,case$kx,case$ky,case$degree)",
    "  bw <- make_bw(case$x,y,case$kx,case$ky,case$degree)",
    "  for (tree in c(FALSE,TRUE)) {",
    "    options(np.tree=tree)",
    "    observed <- npRmpi:::.npcdensbw_eval_only(case$x,data.frame(y=y),bw)$objective[[1L]]",
    "    if (!isTRUE(all.equal(observed,expected,tolerance=2e-10))) stop(sprintf('%s tree=%s observed=%.17g expected=%.17g',case_name,tree,observed,expected))",
    "  }",
    "}",
    "xdup <- data.frame(x=c(0,0,0,1,2)); ydup <- data.frame(y=c(0,0,0,0.5,1))",
    "bwdup <- npcdensbw(xdat=xdup,ydat=ydup,bws=c(2L,2L),bandwidth.compute=FALSE,bwmethod='cv.ml',bwtype='generalized_nn',bwscaling=FALSE)",
    "invalid <- npRmpi:::.npcdensbw_eval_only(xdup,ydup,bwdup,invalid.penalty='dbmax')$objective[[1L]]",
    "stopifnot(identical(invalid,-.Machine$double.xmax))",
    sprintf("cat('%s\\n')", ok_tag)
  )
  result <- npRmpi_run_rscript_subprocess(lines=lines, timeout=90L, env=env)
  info <- paste(result$output, collapse="\n")
  expect_equal(result$status, 0L, info=info)
  expect_true(any(grepl(ok_tag, result$output, fixed=TRUE)), info=info)
})
