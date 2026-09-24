test_that("new fractional cumulative cutoffs use retained support consistently", {
  old <- options(np.messages=FALSE,np.tree=FALSE,np.largeh=FALSE,np.largelambda=FALSE)
  on.exit(options(old),add=TRUE)
  for (s in list(0:2,.212+0:2,c(0,.5,1.7),c(0,2,5))) {
    cuts <- sort(unique(c(min(s)-.25,s,s[-length(s)]+.25,max(s)+.25)))
    ids <- c(1L,2L,1L,2L) # final declared level unused
    d <- data.frame(o=ordered(s[ids],levels=s))
    e <- data.frame(o=ordered(cuts,levels=cuts))
    s <- as.numeric(levels(d$o));cuts <- as.numeric(levels(e$o))
    oracle <- function(lambda,kernel) {
      lattice <- all(abs((s-s[1])-round(s-s[1]))<1e-14)
      support <- if(kernel=="liracine" && lattice)seq(min(s),max(s),by=1)else s
      w <- lambda^abs(outer(s[ids],support,"-"))
      if(kernel=="racineliyan")w <- w/rowSums(w)
      vapply(cuts,function(q)rowSums(w[,support<=q,drop=FALSE]),numeric(length(ids)))
    }
    for(kernel in c("liracine","racineliyan")) {
      for(lambda in c(0,.35,1)) {
        z <- npksum(txdat=d,exdat=e,bws=lambda,okertype=kernel,
          operator="integral",return.kernel.weights=TRUE)
        expect_equal(unname(z$kw),oracle(lambda,kernel),tolerance=2e-12)
        expect_equal(as.numeric(z$ksum),colSums(oracle(lambda,kernel)),tolerance=2e-12)
        expect_true(all(apply(z$kw,1,function(row)all(diff(row)>=-2e-14))))
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

test_that("new off-lattice donors keep the same LR integration interval", {
  old <- options(np.messages=FALSE,np.tree=FALSE,np.largeh=FALSE,np.largelambda=FALSE)
  on.exit(options(old),add=TRUE)
  base <- data.frame(o=ordered(c(0,3),levels=c(0,3)))
  ns <- asNamespace("npRmpi")
  b <- get("kbandwidth.numeric",ns)(bw=.35,xdati=get("untangle",ns)(base),okertype="liracine")
  t <- c(-1.25,.5,1.5,2.25,4.2);q<-c(-.25,0,.25,1,1.25,2,2.25,3,3.25)
  d <- data.frame(o=ordered(t,levels=t));e<-data.frame(o=ordered(q,levels=q))
  oracle <- function(lambda)vapply(q,function(cut)
    rowSums(lambda^abs(outer(t,(0:3)[0:3<=cut],"-"))),numeric(length(t)))
  for(lambda in c(0,.35,1)) {
    b$bw[] <- lambda
    z <- npksum(txdat=d,exdat=e,bws=b,operator="integral",return.kernel.weights=TRUE)
    expect_equal(unname(z$kw),oracle(lambda),tolerance=2e-12)
    expect_true(all(apply(z$kw,1,function(row)all(diff(row)>=-2e-14))))
  }
  b$bw[] <- .35;h<-1e-6
  z<-npksum(txdat=d,exdat=e,bws=b,operator="integral",compute.score=TRUE,
    return.kernel.weights=TRUE,return.derivative.kernel.weights=TRUE)
  expect_equal(as.numeric(z$p.kw),as.numeric((oracle(.35+h)-oracle(.35-h))/(2*h)),tolerance=2e-9)
  # No retained gap is below one here; latent unit-interval points still
  # make the new 1.5 donor's endpoint score singular.
  b$bw[] <- 0
  donor<-data.frame(o=ordered(1.5))
  expect_error(npksum(txdat=donor,exdat=base,bws=b,operator="integral",compute.score=TRUE),
    "score at lambda = 0 is not finite")
  # New cutoffs on large gapped support must not allocate or loop over span.
  large<-data.frame(o=ordered(c(0,1048576),levels=c(0,1048576)))
  query<-data.frame(o=ordered(1048575.25))
  z<-npksum(txdat=large,exdat=query,bws=.5,okertype="liracine",operator="integral",
    return.kernel.weights=TRUE)
  expect_equal(as.numeric(z$kw),c(2,1),tolerance=2e-12)
})
