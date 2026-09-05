test_that("regression gradients retain rows when there is one evaluation point", {
  expect_true(spawn_mpi_slaves(1L))
  on.exit(close_mpi_slaves(force=TRUE),add=TRUE)
  old <- options(np.messages=FALSE)
  on.exit(options(old),add=TRUE)
  set.seed(517)
  n <- 16L
  x <- data.frame(z=runif(n),o=ordered(rep(c("a","b"),8)),u=runif(n))
  y <- x$z+x$u+sin(seq_len(n))
  for(p in 1:3) for(type in c("fixed","generalized_nn","adaptive_nn")) {
    tx <- x[,seq_len(p),drop=FALSE]
    h <- vapply(tx,function(v) if(is.factor(v)) .2 else
      if(type=="fixed") .4 else 8,numeric(1))
    bw <- npregbw(xdat=tx,ydat=y,bws=h,bwtype=type,
      bandwidth.compute=FALSE,regtype="lc")
    one <- npreg(bws=bw,txdat=tx,tydat=y,exdat=tx[1,,drop=FALSE],
      gradients=TRUE,se=TRUE)
    batch <- npreg(bws=bw,txdat=tx,tydat=y,exdat=tx[1:2,,drop=FALSE],
      gradients=TRUE,se=TRUE)
    expect_identical(dim(one$grad),c(1L,p))
    expect_identical(dim(one$gerr),c(1L,p))
    expect_identical(unname(one$grad),unname(batch$grad[1,,drop=FALSE]))
    expect_identical(unname(one$gerr),unname(batch$gerr[1,,drop=FALSE]))
  }
})
