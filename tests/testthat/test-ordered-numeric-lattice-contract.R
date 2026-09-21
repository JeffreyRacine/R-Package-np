test_that("numeric ordered levels require a safely represented integer lattice", {
  old <- options(np.messages=FALSE); on.exit(options(old))
  for(kernel in c("wangvanryzin","liracine","racineliyan")) {
    invalid <- ordered(rep(c(.5,1,1.5),4))
    expect_error(npudensbw(dat=invalid,bws=.3,bandwidth.compute=FALSE,
                          okertype=kernel),"integer")
    valid <- data.frame(o=ordered(rep(c(.5,2.5,6.5),4)))
    expect_s3_class(npudensbw(dat=valid,bws=.3,bandwidth.compute=FALSE,
                             okertype=kernel),"bandwidth")
  }
  for(levels in list(c("1","Inf"),c("1","3000000000"),c("2","1"))) {
    dat <- ordered(rep(levels,4),levels=levels)
    expect_error(npudensbw(dat=dat,bws=.3,bandwidth.compute=FALSE),
                 "finite|increasing|integer|range")
  }
  rank <- ordered(rep(c("low","medium","high"),4),levels=c("low","medium","high"))
  expect_s3_class(npudensbw(dat=rank,bws=.3,bandwidth.compute=FALSE),"bandwidth")
})
