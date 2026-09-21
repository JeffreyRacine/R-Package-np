test_that("auto-off locality starts before numerical preparation", {
  skip_on_cran()
  result <- npRmpi_run_isolated_contract(c(
    "library(npRmpi)",
    "npRmpi.init(nslaves=1L,quiet=TRUE)",
    "options(np.messages=FALSE)",
    "x <- data.frame(x=seq(-1,1,length.out=24)); y <- sin(2*x$x)",
    "b <- npregbw(xdat=x,ydat=y,bws=.4,bandwidth.compute=FALSE)",
    "a <- npreg(b,txdat=x,tydat=y,gradients=TRUE)",
    "d <- npudensbw(dat=x,bws=.4,bandwidth.compute=FALSE)",
    "da <- npudens(d,tdat=x)",
    "options(npRmpi.autodispatch=FALSE)",
    "s <- new.env(); s$n <- 0L",
    "xx <- function() { s$n <- s$n+1L; x }",
    "v <- npreg(b,txdat=xx(),tydat=y,gradients=TRUE)",
    "stopifnot(s$n==1L,isTRUE(all.equal(fitted(a),fitted(v),tolerance=1e-12)))",
    "stopifnot(isTRUE(all.equal(gradients(a),gradients(v),tolerance=1e-12)))",
    "db <- npudens(d,tdat=x)",
    "stopifnot(isTRUE(all.equal(fitted(da),fitted(db),tolerance=1e-12)))",
    "u <- npudensbw(dat=x,nmulti=1L,itmax=2L)",
    "stopifnot(is.finite(u$fval))",
    "k <- npksum(bws=.4,txdat=x,tydat=y)",
    "stopifnot(length(k$ksum)==24L,all(is.finite(k$ksum)))",
    "err <- tryCatch(npudens(d,tdat=rep(NA_real_,24)),error=identity)",
    "stopifnot(inherits(err,'error'),!isTRUE(getOption('npRmpi.local.regression.mode')))",
    "stopifnot(identical(.Call('C_np_set_local_regression_mode',FALSE,PACKAGE='npRmpi'),FALSE))",
    "options(npRmpi.autodispatch=TRUE)",
    "v <- npreg(b,txdat=x,tydat=y,gradients=TRUE)",
    "stopifnot(isTRUE(all.equal(fitted(a),fitted(v),tolerance=1e-12)),mpi.comm.size()==2L)",
    "cat('AUTODISPATCH_LOCAL_ENTRY_OK\\n'); flush.console()"
  ), marker="AUTODISPATCH_LOCAL_ENTRY_OK", timeout=60L)
  if (is.null(result)) skip("No installed npRmpi library for subprocess proof")
  expect_identical(result$status,0L,info=paste(result$output,collapse="\n"))
  expect_true(result$witnessed)
})
