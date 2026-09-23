r21_overlap_reference <- function(x,e,k,bwtype,kernel) {
  ht <- if(bwtype=="fixed") rep(k,length(x)) else
    vapply(seq_along(x),function(i) {
      distances <- abs(x-x[i])
      if(bwtype=="adaptive_nn") distances <- distances[-i]
      sort(distances)[k]
    },0.0)
  he <- if(bwtype=="fixed") rep(k,length(e)) else
    vapply(e,function(v)sort(abs(x-v))[k],0.0)
  radius <- if(kernel=="uniform")1 else sqrt(5)
  base <- if(kernel=="uniform")function(z)rep(.5,length(z)) else
    function(z)3/(4*sqrt(5))*(1-z^2/5)
  raw <- outer(seq_along(x),seq_along(e),Vectorize(function(i,j) {
    lo <- max(x[i]-radius*ht[i],e[j]-radius*he[j])
    hi <- min(x[i]+radius*ht[i],e[j]+radius*he[j])
    if(lo>=hi) return(0)
    integrate(function(t)base((t-x[i])/ht[i])*base((t-e[j])/he[j]),
              lo,hi,rel.tol=1e-12)$value
  }))
  list(raw=raw,normalized=raw/outer(ht,he))
}

test_that("compact NN convolution tree sums retain both bandwidth supports", {
  old <- options(np.messages=FALSE,np.largeh=FALSE,np.largelambda=FALSE)
  on.exit(options(old),add=TRUE)
  set.seed(30921)
  x <- data.frame(x=runif(23))
  e <- data.frame(x=c(.11,.33,.66,.89))
  for(bwtype in c("fixed","generalized_nn","adaptive_nn"))
    for(kernel in c("epanechnikov","uniform")) {
      k <- if(bwtype=="fixed").15 else 8
      ref <- r21_overlap_reference(x$x,e$x,k,bwtype,kernel)
      for(tree in list(FALSE,TRUE,"auto")) for(divide in c(FALSE,TRUE)) {
        options(np.tree=tree)
        f <- npksum(txdat=x,exdat=e,bws=k,bwtype=bwtype,ckertype=kernel,
          operator="convolution",bandwidth.divide=divide,return.kernel.weights=TRUE)
        expect_equal(unname(f$kw),ref$raw,tolerance=2e-11)
        expect_equal(as.numeric(f$ksum),
          colSums(if(divide)ref$normalized else ref$raw),tolerance=2e-11)
      }
    }
})
