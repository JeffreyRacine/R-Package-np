.beta_gnn_identity_contributions <- function(dat, query, k, mapped, order,
                                               lower, upper, cdf) {
  n <- nrow(dat); p <- ncol(dat)
  vapply(seq_len(nrow(query)), function(j) {
    out <- rep(1, n)
    for (d in seq_len(p)) {
      available <- if (mapped) setdiff(seq_len(n), j) else seq_len(n)
      # Public external lookup also caps counts at n-1; preserve its contract.
      lookup <- min(k, n-1L)
      h <- sort(abs(dat[available,d]-query[j,d]))[lookup] * k/lookup
      width <- upper-lower
      z <- (dat[,d]-lower)/width; t <- (query[j,d]-lower)/width
      value <- rep(0, n)
      for (s in seq_len(order/2L)) {
        coefficient <- (-1)^(s+1L)*choose(order/2L,s)
        tau <- (width/h)^2/s
        value <- value + coefficient * if (cdf)
          pbeta(t,1+z*tau,1+(1-z)*tau) else
          dbeta(z,1+t*tau,1+(1-t)*tau)/width
      }
      out <- out*value
    }
    out
  }, numeric(n))
}

test_that("beta GNN full fits preserve training occurrence identity and SEs", {
  old <- options(np.messages=FALSE,np.largeh=FALSE,np.extendednn=TRUE)
  on.exit(options(old),add=TRUE)
  base <- data.frame(x=c(0,.06,.06,.21,.34,.51,.69,.83,1),
                     z=c(.11,.91,.34,.34,.03,.61,.78,1,0))
  for (p in 1:2) for (order in c(2L,4L)) for (k in c(3L,11L))
    for (width in c(1,7)) for (cdf in c(FALSE,TRUE)) {
      dat <- -2+width*base[,seq_len(p),drop=FALSE]
      b <- do.call(if(cdf) npudistbw else npudensbw,
        list(dat=dat,bws=rep(k,p),bandwidth.compute=FALSE,
          bwtype="generalized_nn",ckertype="beta",ckerorder=order,
          ckerbound="fixed",ckerlb=rep(-2,p),ckerub=rep(-2+width,p)))
      for (external in c(FALSE,TRUE)) {
        query <- if(external) dat[c(2L,4L,9L),,drop=FALSE] else dat
        kw <- .beta_gnn_identity_contributions(dat,query,k,!external,order,
                                               -2,-2+width,cdf)
        expected <- colMeans(kw)
        expected.se <- sqrt(colSums(sweep(kw,2L,expected,"-")^2)/
                              (nrow(dat)*(nrow(dat)-1L)))
        for (tree in c(FALSE,TRUE)) {
          options(np.tree=tree)
          args <- list(bws=b,tdat=dat,se=TRUE)
          if(external) args$edat <- query
          fit <- do.call(if(cdf) npudist else npudens,args)
          expect_equal(as.numeric(fitted(fit)),expected,tolerance=4e-12)
          expect_equal(as.numeric(se(fit)),expected.se,tolerance=4e-12)
        }
      }
    }
})

test_that("beta GNN refits and explicit prediction keep distinct query roles", {
  old <- options(np.messages=FALSE,np.tree=FALSE)
  on.exit(options(old),add=TRUE)
  dat <- data.frame(x=c(.02,.07,.16,.31,.47,.64,.82,.96))
  for(cdf in c(FALSE,TRUE)) {
    constructor <- if(cdf) npudistbw else npudensbw
    estimator <- if(cdf) npudist else npudens
    b <- constructor(~x,data=dat,bws=3,bandwidth.compute=FALSE,
      bwtype="generalized_nn",ckertype="beta",ckerbound="fixed",ckerlb=0,ckerub=1)
    fit <- estimator(b)
    train <- colMeans(.beta_gnn_identity_contributions(dat,dat,3,TRUE,2,0,1,cdf))
    external <- colMeans(.beta_gnn_identity_contributions(dat,dat,3,FALSE,2,0,1,cdf))
    expect_equal(as.numeric(fitted(fit)),train,tolerance=4e-12)
    expect_equal(as.numeric(predict(fit,newdata=dat)),external,tolerance=4e-12)
    expect_gt(max(abs(train-external)),1e-5)
  }
})
