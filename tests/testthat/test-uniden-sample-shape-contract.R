test_that("univariate uncertainty uses training sample size", {
  old <- options(np.messages = FALSE); on.exit(options(old))
  X <- seq(.02,.98,length.out = 31)
  Y <- c(.15,.43,.8)
  for (fun in list(npuniden.boundary, npuniden.reflect)) {
    a <- fun(X, Y, h = .18, a = 0, b = 1)
    b <- fun(X, rep(Y,2), h = .18, a = 0, b = 1)
    expect_identical(a$f, b$f[seq_along(Y)])
    expect_equal(a$sd.f, b$sd.f[seq_along(Y)], tolerance = 0)
    expect_equal(a$sd.F, sqrt(abs(a$F*(1-a$F)/length(X))), tolerance = 0)
    expect_equal(b$sd.F, sqrt(abs(b$F*(1-b$F)/length(X))), tolerance = 0)
  }
  a <- npuniden.sc(X,Y,h=.18,lb=0,ub=10)
  b <- npuniden.sc(X,rep(Y,2),h=.18,lb=0,ub=10)
  expect_true(a$solve.QP && b$solve.QP)
  expect_equal(a$f, b$f[seq_along(Y)], tolerance=0)
  expect_equal(a$se.f, b$se.f[seq_along(Y)], tolerance=0)
  expect_equal(a$se.f.sc, b$se.f.sc[seq_along(Y)], tolerance=1e-12)
  expect_equal(a$se.F, sqrt(abs(a$F*(1-a$F)/length(X))), tolerance=0)
  expect_equal(a$se.F.sc, sqrt(abs(a$F.sc*(1-a$F.sc)/length(X))), tolerance=0)
  # Y=X already had the correct denominator in these two owners.
  for (fun in list(npuniden.boundary, npuniden.reflect)) {
    all <- fun(X,h=.18,a=0,b=1)
    some <- fun(X,Y=X[c(2,15,30)],h=.18,a=0,b=1)
    expect_equal(some$sd.f,all$sd.f[c(2,15,30)],tolerance=0)
  }
})
