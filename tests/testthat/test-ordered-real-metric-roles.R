test_that("metric roles retain support and preserve lattice-only family definitions", {
  old <- options(np.messages=FALSE,np.tree=FALSE,np.largeh=FALSE,np.largelambda=FALSE)
  on.exit(options(old),add=TRUE)
  s <- c(0,.5,1.7,pi);ids <- rep(1:3,each=3L)
  d <- data.frame(o=ordered(s[ids],levels=s));y <- seq_along(ids)/length(ids)
  for (kernel in c("liracine","racineliyan")) {
    b <- npregbw(xdat=d,ydat=y,bws=.35,okertype=kernel,bandwidth.compute=FALSE)
    retained <- unserialize(serialize(b,NULL))
    expect_identical(retained$xdati$all.dlev[[1L]],s)
    e <- data.frame(o=ordered(s,levels=s))
    w <- .35^abs(outer(s[ids],s,"-"))
    if (kernel=="racineliyan") w <- w/rowSums(.35^abs(outer(s[ids],s,"-")))
    fit <- npreg(retained,txdat=d,tydat=y,exdat=e)
    expect_equal(as.numeric(fitted(fit)),colSums(w*y)/colSums(w),tolerance=2e-12)
    # A declared unused level remains informative after deleting observations.
    fit <- npreg(retained,txdat=d[-1,,drop=FALSE],tydat=y[-1],exdat=e)
    expect_equal(as.numeric(fitted(fit)),colSums(w[-1,,drop=FALSE]*y[-1])/
                   colSums(w[-1,,drop=FALSE]),tolerance=2e-12)
  }
  for (construct in list(
    function() npregbw(xdat=d,ydat=y,bws=.35,okertype="wangvanryzin",bandwidth.compute=FALSE),
    function() npudensbw(dat=d,bws=.35,bandwidth.compute=FALSE),
    function() npudistbw(dat=d,bws=.35,bandwidth.compute=FALSE),
    function() npcdensbw(xdat=data.frame(y),ydat=d,bws=c(.3,.35),bandwidth.compute=FALSE),
    function() npcdistbw(xdat=data.frame(y),ydat=d,bws=c(.3,.35),bandwidth.compute=FALSE))) {
    expect_error(construct(),"integer distances")
  }
  for (construct in list(
    function() npudensbw(dat=d,bws=.35,okertype="racineliyan",bandwidth.compute=FALSE),
    function() npudistbw(dat=d,bws=.35,okertype="racineliyan",bandwidth.compute=FALSE),
    function() npcdensbw(xdat=data.frame(y),ydat=d,bws=c(.3,.35),oykertype="racineliyan",bandwidth.compute=FALSE),
    function() npcdistbw(xdat=data.frame(y),ydat=d,bws=c(.3,.35),oykertype="racineliyan",bandwidth.compute=FALSE))) {
    expect_no_error(construct())
  }
})

test_that("alphanumeric order and numeric evaluation distances are not rank recoded", {
  old <- options(np.messages=FALSE,np.tree=FALSE)
  on.exit(options(old),add=TRUE)
  d <- data.frame(o=ordered(c("a","b","c","d","e")))
  z <- npksum(txdat=d,bws=.5,return.kernel.weights=TRUE)$kw
  expect_equal(z[1,5],.5^4,tolerance=0)
  expect_equal(z[3,4],.5,tolerance=0)
  d <- data.frame(o=ordered(c(0,1,2),levels=c(0,1,2)))
  e <- data.frame(o=ordered(c(.25,1.25),levels=c(.25,1.25)))
  for (kernel in c("liracine","racineliyan")) {
    z <- npksum(txdat=d,exdat=e,bws=.5,okertype=kernel,return.kernel.weights=TRUE)$kw
    w <- .5^abs(outer(c(0,1,2),c(.25,1.25),"-"))
    if(kernel=="racineliyan")w <- w/rowSums(.5^abs(outer(0:2,0:2,"-")))
    expect_equal(unname(z),w,tolerance=2e-13)
  }
  for (kernel in c("wangvanryzin","nliracine"))
    expect_error(npksum(txdat=d,exdat=e,bws=.5,okertype=kernel),"integer distances")
})
