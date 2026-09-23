test_that("all compact convolution orders preserve mixed operator tree sums", {
  old <- options(np.messages=FALSE,np.largeh=FALSE,np.largelambda=FALSE)
  on.exit(options(old),add=TRUE)
  set.seed(166)
  x <- data.frame(a=runif(48),b=runif(48),
    o=ordered(rep(c(0,1,4),16),levels=c(0,1,4)))
  e <- x[c(2,9,15,21,26,33,39,47),,drop=FALSE]
  Y <- cbind(sin(x$a),x$b)
  W <- cbind(1,seq_len(48)/49)
  for(bwtype in c("fixed","generalized_nn","adaptive_nn"))
    for(order in c(2L,4L,6L,8L)) {
      call <- list(txdat=x,exdat=e,bws=c(if(bwtype=="fixed")c(.09,.13) else c(6,9),.3),
        bwtype=bwtype,ckertype="epanechnikov",ckerorder=order,okertype="racineliyan",
        operator=c("convolution","normal","convolution"),compute.score=TRUE,
        weights=W,tydat=Y,return.kernel.weights=TRUE,return.derivative.kernel.weights=TRUE)
      options(np.tree=FALSE)
      dense <- do.call(npksum,call)
      for(tree in list(TRUE,"auto")) {
        options(np.tree=tree)
        value <- do.call(npksum,call)
        for(field in c("ksum","kw","p.ksum","p.kw"))
          expect_equal(value[[field]],dense[[field]],tolerance=2e-11,
                       info=paste(bwtype,order,tree,field))
      }
    }
})
test_that("primary tree exports equal coordinate products including LOO diagonals", {
  old <- options(np.messages=FALSE,np.largeh=FALSE,np.largelambda=FALSE)
  on.exit(options(old),add=TRUE)
  set.seed(166)
  x <- data.frame(a=runif(48),b=runif(48))
  # Internal divided-row flags belong to the bandwidth-object/default route;
  # numeric dispatch forwards its dots to the bandwidth constructor instead.
  ns <- environment(npksum)
  kbw <- function(cols,h) get("kbandwidth",ns)(bw=h,
    xdati=get("untangle",ns)(x[cols]),xnames=cols,nobs=nrow(x),
    ckertype="epanechnikov")
  single <- list(a=kbw("a",.09),b=kbw("b",.13))
  for(companion in c("normal","integral"))
    for(reverse in c(FALSE,TRUE)) for(normalized in c(FALSE,TRUE))
    for(loo in c(FALSE,TRUE)) {
      e <- if(loo)x else x[c(2,9,15,21,26,33,39,47),,drop=FALSE]
      options(np.tree=FALSE)
      component <- function(col,h,op) do.call(npksum,list(
        txdat=x[col],exdat=e[col],bws=single[[col]],operator=op,
        bandwidth.divide=normalized,.np.internal.bandwidth.divide.weights=normalized,
        return.kernel.weights=TRUE))$kw
      expected <- component("a",.09,"convolution")*component("b",.13,companion)
      ix <- if(reverse)2:1 else 1:2
      args <- list(txdat=x[ix],bws=kbw(names(x)[ix],c(.09,.13)[ix]),
        operator=c("convolution",companion)[ix],leave.one.out=loo,
        bandwidth.divide=normalized,.np.internal.bandwidth.divide.weights=normalized,
        return.kernel.weights=TRUE)
      if(!loo)args$exdat <- e[ix]
      options(np.tree=TRUE)
      f <- do.call(npksum,args)
      expect_equal(unname(f$kw),unname(expected),tolerance=2e-12)
      if(loo)diag(expected) <- 0
      # With bandwidth.divide=FALSE an integral sum includes its h
      # Jacobian; raw exported CDF kernel weights do not include that h.
      factor <- if(!normalized && companion=="integral").13 else 1
      expect_equal(as.numeric(f$ksum),factor*colSums(expected),tolerance=2e-12)
    }
})
