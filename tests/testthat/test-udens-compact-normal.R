test_that("compact normal bands retain the corrected center and complete variance", {
  x <- qbeta((1:20 - .35)/20.3, 1.4, 2.1)
  grid <- c(.04, .3, .7, .96, .3)
  for (kernel in c("uniform", "epanechnikov")) {
    z <- .np_udens_compact_normal_moments(x, .16, grid, kernel)
    old <- .np_udens_compact_prepare(x, .16, grid, kernel, FALSE)
    expect_identical(z$center, old$center)
    expect_identical(z$se[2L], z$se[5L])
    expect_true(all(is.finite(z$se) & z$se > 0))
    expect_gt(max(abs(z$theta - z$thetaQ)), 1e-5)
    # Independently contract the complete influence on the split rule.
    plan <- .np_udens_compact_operator_plan(grid, .16, kernel)
    r <- .np_udens_compact_rule(z$theta, 128L, plan)
    s <- .np_udens_compact_rule(z$thetaQ, 128L, plan)
    p <- .np_udens_q1_ratio_weights(x, z$theta, r$A, r)$prob
    A <- .np_udens_compact_operator(plan, exp(r$lv))
    Aq <- .np_udens_compact_operator(plan, exp(s$lv))
    density <- exp(drop(cbind(log(grid), log1p(-grid)) %*% z$thetaQ) - s$A)
    G <- sweep(sweep(cbind(log(grid), log1p(-grid)), 2L, s$mean), 1L, density, "*")
    D <- G - Aq %*% (sweep(s$T, 2L, s$mean) * s$prob)
    raw <- A + (D %*% solve(s$I)) %*% t(r$T)
    oracle <- vapply(seq_along(grid), function(j) {
      mu <- sum(p * raw[j,]); sum(p * (raw[j,] - mu)^2)/length(x)
    }, 0)
    expect_equal(z$se^2, oracle, tolerance=1e-11)
    only.kernel <- rowSums((A - drop(A %*% p))^2 * rep(p, each=nrow(A)))/length(x)
    expect_gt(max(abs(oracle - only.kernel)), 1e-5)
    expect_lte(z$workspace, 192*1024^2)
  }
  expect_error(.np_udens_compact_normal_moments(x,.16,0,"uniform"),"invalid")
  expect_error(.np_udens_compact_normal_moments(x,.16,.5,"gaussian"),"invalid")
  expect_error(.np_udens_compact_normal_moments(rep(.5,10000),.16,rep(.5,1000),"uniform"),"workspace")
})

test_that("compact normal public assembly retains units and ordinary critical values", {
  x <- seq(.05,.95,length.out=12L); grid <- c(.2,.7,.2)
  for (kernel in c("uniform","epanechnikov")) {
    z <- .np_udens_compact_normal_moments(x,.2,grid,kernel)
    state <- list(q=1L,n=length(x),variant="asymptotic",normal.pmzsd=TRUE,
      x=data.frame(x=3+4*x),normalized=matrix(x,ncol=1L),h=.2,physical.h=.8,
      lower=3,width=4,jacobian=4,kernel=kernel)
    for (alpha in c(.1,.05,.01)) {
      out <- .np_udens_compact_public_errors(state,data.frame(x=3+4*grid),alpha,7L,"pmzsd")
      expect_equal(out$inference$center,z$center/4,tolerance=1e-12)
      expect_equal(out$inference$se,z$se/4,tolerance=1e-12)
      expect_equal(out$inference$bands$pmzsd[,1L,1L],(z$center-qnorm(1-alpha/2)*z$se)/4,tolerance=1e-12)
      expect_equal(out$inference$bands$pmzsd[,2L,1L],(z$center+qnorm(1-alpha/2)*z$se)/4,tolerance=1e-12)
      expect_identical(out$inference$settings$B,0L)
      expect_false(out$inference$normal.band$calibrated.envelope)
    }
    expect_error(.np_udens_compact_public_errors(state,data.frame(x=3),.05,7L,"pmzsd"),"interior")
  }
})
