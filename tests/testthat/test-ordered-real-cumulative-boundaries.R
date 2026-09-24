test_that("new fractional cumulative cutoffs use retained support consistently", {
  old <- options(np.messages=FALSE,np.tree=FALSE,np.largeh=FALSE,np.largelambda=FALSE)
  on.exit(options(old),add=TRUE)
  for (s in list(0:2,.212+0:2,c(0,.5,1.7),c(0,2,5))) {
    cuts <- sort(unique(c(min(s)-.25,s[-length(s)]+.25,max(s)+.25)))
    ids <- c(1L,2L,1L,2L) # final declared level unused
    d <- data.frame(o=ordered(s[ids],levels=s))
    e <- data.frame(o=ordered(cuts,levels=cuts))
    s <- as.numeric(levels(d$o));cuts <- as.numeric(levels(e$o))
    oracle <- function(lambda,kernel) {
      w <- lambda^abs(outer(s[ids],s,"-"))
      if(kernel=="racineliyan")w <- w/rowSums(w)
      vapply(cuts,function(q)rowSums(w[,s<=q,drop=FALSE]),numeric(length(ids)))
    }
    for(kernel in c("liracine","racineliyan")) {
      for(lambda in c(0,.35,1)) {
        z <- npksum(txdat=d,exdat=e,bws=lambda,okertype=kernel,
          operator="integral",return.kernel.weights=TRUE)
        expect_equal(unname(z$kw),oracle(lambda,kernel),tolerance=2e-12)
        expect_equal(as.numeric(z$ksum),colSums(oracle(lambda,kernel)),tolerance=2e-12)
      }
      h <- 1e-6
      z <- npksum(txdat=d,exdat=e,bws=.35,okertype=kernel,
        operator="integral",compute.score=TRUE,return.kernel.weights=TRUE,
        return.derivative.kernel.weights=TRUE)
      expect_equal(as.numeric(z$p.kw),
        as.numeric((oracle(.35+h,kernel)-oracle(.35-h,kernel))/(2*h)),tolerance=2e-9)
    }
  }
})
