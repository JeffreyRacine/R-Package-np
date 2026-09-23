test_that("partially linear hat fanout retains response names and output shape", {
  skip_if_not(.mpi_pool_active(), "requires a live MPI pool")
  old <- options(np.messages=FALSE, np.tree=FALSE,
    npRmpi.hat.operator.fanout=TRUE, npRmpi.hat.operator.fanout.min.work=0)
  on.exit(options(old), add=TRUE)
  set.seed(827)
  n <- 32L
  z <- data.frame(z=sort(runif(n)))
  x <- data.frame(x=rnorm(n))
  y <- sin(4*z$z)+x$x
  bw <- npplregbw(xdat=x, zdat=z, ydat=y, bws=matrix(.3,2L,1L),
    regtype="lp", degree=2L, bandwidth.compute=FALSE)
  for (training in c(FALSE, TRUE)) {
    ex <- if (training) x else data.frame(x=seq(-.4,.4,length.out=12L))
    ez <- if (training) z else data.frame(z=seq(.2,.8,length.out=12L))
    args <- list(bws=bw, txdat=x, tzdat=z, exdat=ex, ezdat=ez, output="apply")
    for (labels in list(NULL, c("response","shifted"), c("same","same"), c("","second"))) {
      yy <- cbind(y, y+1)
      colnames(yy) <- labels
      args$y <- yy
      options(npRmpi.hat.operator.fanout=FALSE)
      local <- do.call(npplreghat, args)
      options(npRmpi.hat.operator.fanout=TRUE)
      parallel <- do.call(npplreghat, args)
      expect_identical(colnames(parallel), colnames(yy))
      expect_identical(dim(parallel), c(nrow(ex), 2L))
      expect_equal(parallel, local, tolerance=3e-12)
    }
    for (yy in list(y, matrix(y,ncol=1L,dimnames=list(NULL,"one")))) {
      args$y <- yy
      options(npRmpi.hat.operator.fanout=FALSE)
      local <- do.call(npplreghat, args)
      options(npRmpi.hat.operator.fanout=TRUE)
      parallel <- do.call(npplreghat, args)
      expect_null(dim(parallel))
      expect_null(names(parallel))
      expect_equal(parallel, local, tolerance=3e-12)
    }
    args$y <- cbind(response=y,shifted=y+1)
    applied <- do.call(npplreghat,args)
    args$output <- "matrix"
    H <- do.call(npplreghat,args)
    expect_equal(applied,H %*% args$y,tolerance=3e-12)
  }
})
