test_that("MPI generic LC hats retain the deleted neighbor geometry", {
  old <- options(np.messages=FALSE,np.extendednn=TRUE,np.tree=FALSE)
  on.exit(options(old),add=TRUE)
  x <- c(.03,.09,.18,.29,.43,.57,.71,.84,.97)
  d <- data.frame(x=x);y <- sin(9*x)+x;n <- length(x)
  for(type in c("generalized_nn","adaptive_nn")) for(kernel in c("gaussian","epanechnikov")) {
    k <- n+2L;kk <- n-2L
    w <- vapply(seq_len(n),function(i) {
      donor <- setdiff(seq_len(n),i)
      h <- if(type=="generalized_nn")sort(abs(x[donor]-x[i]))[kk]*k/kk else
        vapply(donor,function(j)sort(abs(x[setdiff(donor,j)]-x[j]))[kk]*k/kk,0.)
      u <- (x[donor]-x[i])/h
      out <- numeric(n)
      out[donor] <- (if(kernel=="gaussian")dnorm(u) else
        .75*(1-u^2/5)*(abs(u)<=sqrt(5))/sqrt(5))/h
      out
    },numeric(n))
    expected <- t(sweep(w,2L,colSums(w),"/"))
    b <- npregbw(xdat=d,ydat=y,bws=k,bwtype=type,regtype="lc",
      ckertype=kernel,bandwidth.compute=FALSE)
    H <- npreghat(b,txdat=d,leave.one.out=TRUE)
    expect_equal(unname(as.matrix(H)),unname(expected),ignore_attr=TRUE,tolerance=2e-13)
    expect_equal(as.numeric(npreghat(b,txdat=d,y=y,output="apply",leave.one.out=TRUE)),
      as.vector(expected%*%y),tolerance=2e-13)
    expect_equal(unname(npreghat(b,txdat=d,y=y,output="constraint",leave.one.out=TRUE)),
      unname(t(expected)*y),ignore_attr=TRUE,tolerance=2e-13)
    raw <- npksum(txdat=d,bws=b,leave.one.out=TRUE,return.kernel.weights=TRUE,
      bandwidth.divide=TRUE)$kw
    expect_true(max(abs(t(sweep(raw,2L,colSums(raw),"/"))-expected))>1e-6)
  }
})
