test_that("fixed beta CDF bands condition on the declared continuous support", {
  x <- seq(2, 5, length.out = 45)^1; x <- x[-c(1,45)]
  dat <- data.frame(x)
  b <- npudistbw(dat=dat,bws=.47,bandwidth.compute=FALSE,
                 ckertype="beta",ckerbound="fixed",ckerlb=2,ckerub=5)
  st <- .np_udist_bc_prepare(b,dat,allow.compact=TRUE)
  g <- c(5,2.6,3.2,2,4.5)
  evaluate <- function(g, bands=c("pointwise","simultaneous","bonferroni"))
    .np_udist_joint_inference(st,data.frame(x=g),B=199L,bands=bands)
  set.seed(12); z <- evaluate(g)
  inner <- c(2,3,5); set.seed(12); zi <- evaluate(g[inner])
  for(name in names(z$intervals)){
    expect_equal(z$intervals[[name]][inner,,1],zi$intervals[[name]][,,1],tolerance=1e-13)
    expect_equal(unname(z$intervals[[name]][c(1,4),,1]),cbind(c(1,0),c(1,0)))
  }
  h <- .47/3; u <- (x-2)/3
  component <- function(k) vapply((g-2)/3,function(t)
    pbeta(t,1+u/(k*h^2),1+(1-u)/(k*h^2)),numeric(length(x)))
  J <- 2*component(1L)-component(2L)
  expect_equal(z$center,colMeans(J),tolerance=2e-13)
  expect_equal(z$se,sqrt(colSums(sweep(J,2,colMeans(J))^2))/length(x),tolerance=2e-13)
  expect_identical(z$settings$structural.endpoints,c(1L,4L))
  expect_equal(z$settings$stochastic.family.size,3)
  expect_equal(z$raw.max,zi$raw.max,tolerance=1e-13)
  entry <- .Random.seed; endpoint <- evaluate(c(2,5,2))
  expect_identical(.Random.seed,entry)
  expect_null(endpoint$raw.max)
  expect_identical(endpoint$settings$draws,0L)
  expect_equal(endpoint$se,c(0,0,0))
  set.seed(13); single <- evaluate(c(2,3,5))
  expect_identical(single$intervals$simultaneous,single$intervals$pointwise)
  expect_identical(single$intervals$bonferroni,single$intervals$pointwise)
  expect_error(evaluate(5.1),"within the declared support")
  tight <- .np_udist_beta_support(b,x,.01)
  st$beta<-tight
  expect_error(evaluate(5-1e-12),"interior variation")
})

