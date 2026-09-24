test_that("raw folded NN geometry is independent of preceding scaling state", {
  skip_on_cran()
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force=TRUE),add=TRUE)
  old <- options(np.messages=FALSE,np.extendednn=TRUE,np.tree=FALSE)
  on.exit(options(old),add=TRUE)
  x <- data.frame(x=c(.02,.05,.12,.18,.4,.56,.68,.87,.98),
    u=factor(rep(c("a","b","c"),3)),o=ordered(rep(1:3,3)))
  y <- sin(7*x$x);n <- nrow(x)
  # First native operation must be beta ANN in a fresh qualification process.
  # The original raw-fold owner consulted a NULL ambient standard deviation.
  for(kernel in c("beta","gaussian"))for(type in c("adaptive_nn","generalized_nn"))
    for(count in c(n+2L,3L)) {
      b <- npregbw(xdat=x,ydat=y,bws=c(count,.2,.3),bwtype=type,
        bandwidth.compute=FALSE,regtype="lc",ckertype=kernel,
        ckerbound="fixed",ckerlb=0,ckerub=1)
      actual <- .np_kernel_weights_direct(b,x,leave.one.out=TRUE)
      expected <- vapply(seq_len(n),function(i) {
        out <- numeric(n)
        literal <- npksum(txdat=x[-i,,drop=FALSE],exdat=x[i,,drop=FALSE],
          tydat=diag(n-1L),bws=b,bandwidth.divide=TRUE,return.kernel.weights=TRUE)
        # Legacy GNN exports raw kernels; query normalization cancels in hats.
        # Beta and ANN internal rows carry their absolute kernel weights.
        out[-i] <- if(kernel!="beta" && type=="generalized_nn")literal$kw else literal$ksum
        out
      },numeric(n))
      expect_equal(unname(actual),expected,ignore_attr=TRUE,tolerance=2e-12,
        info=paste(kernel,type,count))
      # MPI hat folding is the separate H adapter contract, not H2 proof.
      for(scale in c(TRUE,FALSE)) {
        warm <- npregbw(xdat=x$x,ydat=y,bws=.4,bwscaling=scale,bandwidth.compute=FALSE)
        invisible(npreg(warm,txdat=x$x,tydat=y))
        expect_identical(.np_kernel_weights_direct(b,x,leave.one.out=TRUE),actual)
      }
      # Public non-fold leave.one.out is the accepted zero-diagonal convention.
      ordinary <- matrix(as.numeric(npksum(txdat=x,tydat=diag(n),bws=b,
        bandwidth.divide=TRUE)$ksum),nrow=n)
      zeroed <- matrix(as.numeric(npksum(txdat=x,tydat=diag(n),bws=b,
        leave.one.out=TRUE,bandwidth.divide=TRUE)$ksum),nrow=n)
      expected.zero <- ordinary;diag(expected.zero)<-0
      expect_equal(zeroed,expected.zero,tolerance=2e-12,ignore_attr=TRUE)
      expect_identical(.np_kernel_weights_direct(b,x,leave.one.out=TRUE),actual)
    }
})
