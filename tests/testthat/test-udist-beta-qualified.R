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


test_that("beta CDF inference refuses boundary observations and encloses signed contributions", {
  dat <- data.frame(x = seq(.04, .96, length.out = 41))
  bw <- npudistbw(dat = dat, bws = .12, bandwidth.compute = FALSE,
                  ckertype = "beta", ckerbound = "fixed", ckerlb = 0, ckerub = 1)
  for (value in c(0, 1, -.01, 1.01)) {
    altered <- dat
    altered$x[1] <- value
    expect_error(.np_udist_bc_prepare(bw, altered, allow.compact = TRUE),
                 "strictly inside|support|bound")
  }
  support <- .np_udist_beta_support(bw, dat$x, .12)
  targets <- c(.06, .3, .65, .94)
  actual <- .np_udist_beta_ranges(targets, support)
  contributions <- function(u, q, scale)
    pbeta(q, 1 + u/(scale*.12^2), 1 + (1-u)/(scale*.12^2))
  cells <- seq(0, 1, length.out = 65)
  dense <- seq(0, 1, length.out = 2049)
  for (i in seq_along(targets)) {
    one <- contributions(cells, targets[i], 1)
    two <- contributions(cells, targets[i], 2)
    expect_equal(actual$lower[i], max(-1, min(2*one[-1]-two[-65])-2e-13), tolerance = 1e-13)
    expect_equal(actual$upper[i], min(2, max(2*one[-65]-two[-1])+2e-13), tolerance = 1e-13)
    j <- 2*contributions(dense, targets[i], 1)-contributions(dense, targets[i], 2)
    expect_true(all(actual$lower[i] <= j & j <= actual$upper[i]))
  }
  st <- .np_udist_bc_prepare(bw, dat, allow.compact = TRUE)
  z <- .np_udist_joint_inference(st, data.frame(x = c(0, 1)), B = 99L)
  for (interval in z$intervals) expect_identical(as.vector(interval), c(0,1,0,1))
})
