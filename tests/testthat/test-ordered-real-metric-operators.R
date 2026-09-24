# Independent finite counting-support definitions, not package kernels.
metric179_oracle <- function(s, lambda, kernel, op) {
  w <- lambda^abs(outer(s,s,"-"))
  if (kernel == "racineliyan") w <- w/rowSums(w)
  switch(op,normal=w,integral=t(apply(w,1L,cumsum)),convolution=tcrossprod(w))
}

test_that("real ordered LR and RLY operators use the original metric and retained support", {
  old <- options(np.messages=FALSE,np.tree=FALSE,np.largeh=FALSE,np.largelambda=FALSE)
  on.exit(options(old),add=TRUE)
  supports <- list(c(0,.5,1.7,pi), .212+c(0,.5,1.7,pi),
                   c(0,1/1001,2), c(-2.5,0,1.25,1e12),
                   c(0,.Machine$integer.max+10))
  for (s in supports) {
    ids <- rep(seq_len(length(s)-1L),2L) # final declared level is unused
    d <- data.frame(o=ordered(s[ids],levels=s))
    e <- data.frame(o=ordered(s,levels=s))
    for (kernel in c("liracine","racineliyan"))
      for (lambda in c(0,.35,1)) for (op in c("normal","integral","convolution")) {
        expected <- metric179_oracle(s,lambda,kernel,op)[ids,,drop=FALSE]
        z <- npksum(txdat=d,exdat=e,bws=lambda,okertype=kernel,operator=op,
                    return.kernel.weights=TRUE)
        expect_equal(unname(z$kw),expected,tolerance=2e-12)
        expect_equal(as.numeric(z$ksum),colSums(expected),tolerance=2e-12)
      }
  }
})

test_that("real ordered scores differentiate normal and paired operators", {
  old <- options(np.messages=FALSE,np.tree=FALSE,np.largeh=FALSE,np.largelambda=FALSE)
  on.exit(options(old),add=TRUE)
  s <- c(0,.5,1.7,pi);ids <- c(1,2,1,3,2,3)
  d <- data.frame(o=ordered(s[ids],levels=s));e <- data.frame(o=ordered(s,levels=s))
  h <- 1e-6
  for (kernel in c("liracine","racineliyan"))
    for (op in c("normal","integral","convolution")) {
      expected <- (metric179_oracle(s,.35+h,kernel,op)-
                   metric179_oracle(s,.35-h,kernel,op))/(2*h)
      z <- npksum(txdat=d,exdat=e,bws=.35,okertype=kernel,operator=op,
                  compute.score=TRUE,return.kernel.weights=TRUE,
                  return.derivative.kernel.weights=TRUE)
      expect_equal(as.numeric(z$p.kw),as.numeric(expected[ids,,drop=FALSE]),tolerance=2e-9)
    }
  for (kernel in c("liracine","racineliyan"))
    expect_error(npksum(txdat=d,bws=0,okertype=kernel,compute.score=TRUE),
                 "score at lambda = 0 is not finite")
  # Nonintegral gaps above one have a finite zero score at the origin.
  s <- c(0,1.5,3.25);d <- data.frame(o=ordered(s,levels=s))
  for (kernel in c("liracine","racineliyan")) {
    z <- npksum(txdat=d,bws=0,okertype=kernel,compute.score=TRUE,
                return.kernel.weights=TRUE,return.derivative.kernel.weights=TRUE)
    expect_equal(as.numeric(z$p.kw),rep(0,9),tolerance=0)
  }
})
