test_that("univariate uncertainty uses training sample size", {
  skip_on_cran()
  owns.pool <- !.mpi_pool_active()
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  if (owns.pool) on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
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

test_that("one-sided density constraints build conformable quadratic programs", {
  skip_on_cran()
  owns.pool <- !.mpi_pool_active()
  if (!spawn_mpi_slaves()) skip("Could not spawn MPI slaves")
  if (owns.pool) on.exit(close_mpi_slaves(force = TRUE), add = TRUE)
  X <- seq(.05,.95,length.out=12)
  for (bounds in list(list(lb=0),list(ub=10),list(lb=0,ub=10),list(lb=0,ub=Inf),
                      list(lb=.9),list(ub=1.1))) {
    fit <- do.call(npuniden.sc,c(list(X=X,Y=X,h=.25),bounds))
    expect_true(fit$solve.QP)
    if (!is.null(bounds$lb)) expect_gte(min(fit$f.sc)+1e-10,bounds$lb)
    if (!is.null(bounds$ub)) expect_lte(max(fit$f.sc)-1e-10,bounds$ub)
  }
  expect_error(npuniden.sc(X,h=.25,lb=c(0,.1)),"numeric scalars")
})
