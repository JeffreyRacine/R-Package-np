test_that("bounded normal admission is separate from calibrated bands", {
  b <- list(ndim=1L)
  admit <- .np_udens_bounded_normal_request
  expect_true(admit(b,"asymptotic","pmzsd",NULL,.05))
  for (band in c("pointwise","simultaneous","bonferroni","all"))
    expect_false(admit(b,"asymptotic",band,NULL,.05))
  expect_error(admit(b,"bootstrap","pmzsd",NULL,.05),"asymptotic")
  b$ndim <- 2L
  expect_error(admit(b,"asymptotic","pmzsd",NULL,.05),"q=1")
  b$ndim <- 1L
  expect_error(admit(b,"asymptotic","pmzsd","symmetric-sd",.05),"calibration")
  for (alpha in list(0,.5,NA_real_,Inf,c(.1,.05)))
    expect_error(admit(b,"asymptotic","pmzsd",NULL,alpha),"alpha")
})

test_that("normal reference fitting uses the beta score equation", {
  st <- .np_udens_bounded_normal_analytic(c(0,0))
  expect_equal(st$A,0,tolerance=1e-14)
  expect_equal(st$mean,c(-1,-1),tolerance=1e-14)
  expect_equal(st$I,matrix(c(1,1-pi^2/6,1-pi^2/6,1),2),tolerance=1e-14)
  z <- .np_udens_bounded_normal_fit(c(-1,-1),c(.2,.1))
  expect_lt(max(abs(z$theta)),2e-8)
  expect_lt(z$score_error,2e-8)
  expect_error(.np_udens_bounded_normal_analytic(c(-1,0)),"reference guide")
})

test_that("normal reference stopping satisfies both score certificates", {
  target <- c(-0.9517478719711074, -0.8496103883630108)
  start <- c(0.22875755046493026, 0.43028131723618307)
  premature <- c(0.36439115274809775, 0.4653588189532693)
  score <- function(theta) digamma(theta + 1) - digamma(sum(theta + 1))
  limit <- 1e-12 + 2e-8 * max(abs(target))
  expect_lt(max(abs(score(premature) - target)), 2e-8)
  expect_gt(max(abs(score(premature) - target)), limit)
  z <- .np_udens_bounded_normal_fit(target, start)
  expect_lte(max(abs(score(z$theta) - target)), limit)
  expect_lt(z$score_error, 1e-12)
  for (shape in list(c(.8, 1.2), c(1, 1), c(4, 9), c(50, 50))) {
    goal <- digamma(shape) - digamma(sum(shape))
    z <- .np_udens_bounded_normal_fit(goal, shape * 1.05 - 1)
    residual <- max(abs(score(z$theta) - goal))
    expect_lt(residual, 2e-8)
    expect_lte(residual, 1e-12 + 2e-8 * max(abs(goal)))
    expect_true(all(z$theta > -1))
  }
})

test_that("complete normal variance keeps guide and cross contributions", {
  A <- rbind(c(1,3,2,5),c(2,-1,4,1))
  T <- cbind(c(-2,-1,0,1),c(1,-1,2,0)); p <- c(.1,.2,.3,.4)
  D <- rbind(c(.5,-.3),c(-.2,.4)); I <- matrix(c(2,.3,.3,1),2)
  z <- .np_udens_bounded_normal_variance(A,T,p,D,I,20)
  oracle <- vapply(1:2,function(j) {
    coefficients <- solve(I,D[j,])
    u <- A[j,]+T[,1]*coefficients[1]+T[,2]*coefficients[2]
    sum(p*(u-sum(p*u))^2)/20
  },0)
  expect_equal(z$value,oracle,tolerance=1e-14)
  drop <- .np_udens_bounded_normal_variance(A,T,p,D*0,I,20)
  expect_gt(max(abs(drop$value-oracle)),.001)
  expect_error(.np_udens_bounded_normal_budget(1e8,128,50),"workspace")
})

test_that("bounded normal output preserves center and physical support units", {
  x <- seq(.06,.94,length.out=16L); grid <- c(.3,.7,.3)
  for (kernel in c("beta","gaussian")) {
    z <- .np_udens_bounded_normal_moments(x,.15,grid,kernel)
    score <- if(kernel=="beta") .np_udens_beta_q1_score else .np_udens_fgd_score
    old <- .np_udens_beta_q1_moments(x,.15,grid,z$theta,
      if(kernel=="beta")64L else 24L,128L,FALSE,score)
    expect_equal(z$center,old$center,tolerance=1e-12)
    state <- list(q=1L,n=length(x),variant="asymptotic",x=data.frame(x=3+4*x),
      normalized=matrix(x,ncol=1L),h=.15,lower=3,width=4,jacobian=4)
    out <- .np_udens_bounded_normal_public(state,data.frame(x=3+4*grid),.05,"pmzsd",kernel)
    expect_equal(out$inference$center,z$center/4,tolerance=1e-12)
    expect_equal(out$inference$se,z$se/4,tolerance=1e-12)
    expect_equal(out$inference$bands$pmzsd[,1,1],(z$center-qnorm(.975)*z$se)/4,tolerance=1e-12)
    expect_equal(out$inference$bands$pmzsd[,2,1],(z$center+qnorm(.975)*z$se)/4,tolerance=1e-12)
    expect_false(out$inference$normal.band$calibrated.envelope)
    expect_identical(out$inference$settings$B,0L)
    expect_null(out$inference[["pilot.residual.scale",exact=TRUE]])
    sliced <- .np_udens_beta_plot_subset(out,c(3L,2L))
    expect_identical(sliced$inference$se,out$inference$se[c(3L,2L)])
    expect_error(.np_udens_bounded_normal_public(state,data.frame(x=3),.05,"pmzsd",kernel),"interior")
  }
})
