test_that('q1 variance planning preserves literal target order and guide measure', {
  theta <- c(.043,.238); grid <- c(.0027,.206,.9985)
  p <- .np_udens_q1_variance_plan(theta,.00688,grid)
  expect_identical(p,.np_udens_q1_variance_plan(theta,.00688,c(rev(grid),grid)))
  # Distinct targets in the same logit interval add no new panel boundaries.
  dense <- plogis(seq(-1.9,-1.1,length.out=10001L))
  expect_identical(.np_udens_q1_variance_plan(theta,.00688,dense),
    .np_udens_q1_variance_plan(theta,.00688,plogis(-1.5)))
  z <- .np_udens_q1_variance_rule(theta,p)
  ex <- .np_udens_bounded_normal_analytic(theta)
  expect_true(all(z$prob>0))
  expect_equal(sum(z$prob),1,tolerance=1e-14)
  expect_equal(z$A,ex$A,tolerance=1e-12)
  expect_equal(unname(z$mean),ex$mean,tolerance=1e-12)
  expect_equal(unname(z$I),ex$I,tolerance=1e-11)
  expect_lte(z$nodes,8192L)
  expect_error(.np_udens_q1_variance_plan(theta,.00688,c(0,.5)),'geometry')
  expect_error(.np_udens_q1_variance_plan(theta,1e-8,grid),'node budget')
  p$nodes <- p$nodes-1L
  expect_error(.np_udens_q1_variance_rule(theta,p),'variance plan')
})

test_that('q1 batched columns retain independent score arithmetic', {
  grid <- c(.1,.4,.9,.4); x <- seq(.03,.97,length.out=12)
  coordinates <- list(list(lv=log(x),l1v=log1p(-x)),
    list(lv=log(rev(x)),l1v=log1p(-rev(x))))
  for(score in list(.np_udens_beta_q1_score,.np_udens_fgd_score)) {
    z <- .np_udens_q1_score_columns(grid,coordinates,.08,.2,32L,score)
    for(i in 1:2) expect_identical(z[[i]],
      score(grid,coordinates[[i]]$lv,coordinates[[i]]$l1v,.08*.08,.2,32L))
  }
})

test_that('q1 common positive quadrature keeps complete Gram covariance', {
  x <- seq(.03,.97,length.out=24); g <- c(.21,.5,.79)
  theta <- .np_udens_beta_q1_fit(x)$theta
  a <- .np_udens_beta_q1_moments(x,.08,g,theta,64L,128L,TRUE)
  b <- .np_udens_beta_q1_moments(x,.08,g,theta,128L,256L,TRUE)
  d <- .np_udens_beta_q1_moments(x,.08,g,theta,64L,128L,FALSE)
  expect_identical(a$pilot_variance,d$pilot_variance)
  expect_equal(diag(a$covariance),a$pilot_variance,tolerance=1e-14)
  expect_gte(min(eigen(a$covariance,symmetric=TRUE,only.values=TRUE)$values),-1e-10)
  expect_silent(.np_udens_beta_check(a$covariance,b$covariance,'fixture covariance'))
  expect_error(.np_udens_beta_check(a$covariance,b$covariance+.01,'wrong covariance'),'refinement failed')
  p <- .np_udens_q1_variance_plan(theta,.08,g)
  z <- .np_udens_q1_variance_rule(theta,p)
  expect_error(.np_udens_beta_check(sum(z$prob*1.01),1,'wrong mass'),'refinement failed')
})

test_that('q1 squared-score peaks meet the unchanged refinement certificate', {
  # Compact geometry from the saved narrow-bandwidth failure. The predecessor
  # misses this variance by >9e-8 against an 8.3e-9 limit; no large data fixture.
  x <- seq(.03,.97,length.out=24)
  h <- .0068793233388765741
  grid <- .20589632514546927
  theta <- c(.043018331478746313,.23782979033327017)
  a <- .np_udens_beta_q1_moments(x,h,grid,theta,64L,128L,TRUE)
  b <- .np_udens_beta_q1_moments(x,h,grid,theta,128L,256L,TRUE)
  for(field in c('pilot_variance','guide_variance','covariance'))
    expect_silent(.np_udens_beta_check(a[[field]],b[[field]],field))
})
