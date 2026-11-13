test_that("Q6 compact operators retain their reproducing moments", {
  for (kernel in c("uniform", "epanechnikov")) {
    span <- if (kernel == "uniform") 1 else sqrt(5)
    for (k in 0:5) {
      value <- integrate(function(r) r^k * .np_udens_q6_kernel(r, kernel),
                         -span, span, abs.tol = 1e-12, rel.tol = 1e-11)$value
      expect_equal(value, as.numeric(k == 0), tolerance = 1e-11)
    }
  }
})

test_that("Q6 uses the complete Student6 reference and original scale", {
  x <- qnorm((1:20)/21)
  h <- .6
  grid <- c(-.7, .2, .8)
  for (kernel in c("uniform", "epanechnikov")) {
    covariance <- .np_udens_q6_covariance(x, h, grid, kernel, 96L)
    v <- mean((x-mean(x))^2)
    expect_equal(covariance$m4, mean((x-mean(x))^4) + 6*h*h*v + 6*h^4)
    ref <- .np_udens_q6_reference(x, h, grid, kernel)
    expect_equal(ref$se, sqrt(diag(covariance$C)), tolerance = 1e-12)
    expect_true(all(ref$diagnostics <= 1))
    # Independently integrate the reference density against the compact operator.
    span <- if (kernel == "uniform") 1 else sqrt(5)
    q <- function(u) vapply(u, function(t) mean(dt((t-x)/(h*sqrt(2/3)),6))/(h*sqrt(2/3)), 0.0)
    oracle <- vapply(grid, function(t) {
      integrate(function(r) .np_udens_q6_kernel(r,kernel) *
        (q(t+h*r)-dnorm(t+h*r,mean(x),sqrt(v+h*h))), -span, span,
        rel.tol=1e-11, abs.tol=1e-12)$value + dnorm(t,mean(x),sqrt(v+h*h))
    }, 0.0)
    expect_equal(ref$reference.center, oracle, tolerance = 1e-10)
    expect_false(isTRUE(all.equal(ref$center, ref$reference.center, tolerance=1e-9)))
    values <- c(.5, 1, 2, 3)
    out <- .np_udens_q6_intervals(ref, values, .1)
    critical <- unname(quantile(values,.9,type=7))
    expect_identical(out$critical, critical)
    expect_identical(out$lower, ref$center-critical*ref$se)
    expect_identical(out$upper, ref$center+critical*ref$se)
    perm <- .np_udens_q6_reference(x,h,grid[c(3,1,3,2)],kernel)
    expect_identical(perm$covariance, ref$covariance)
    expect_identical(.np_udens_q6_intervals(perm,values,.1)$lower,
                     out$lower[c(3,1,3,2)])
  }
})

test_that("Q6 numerical failures remain failures", {
  expect_error(.np_udens_q6_geometry(1:9,.1,0,"uniform"), "invalid Q6")
  expect_error(.np_udens_q6_geometry(rep(1,10),.1,0,"uniform"), "guide variance")
  expect_error(.np_udens_q6_geometry(1:10,0,0,"uniform"), "invalid Q6")
  expect_error(.np_udens_q6_budget(1000,100000,9999), "workspace bound")
  expect_error(.np_udens_q6_check(1,1.001,"test"), "refinement failed")
  expect_error(.np_udens_q6_check(Inf,1,"test"), "nonfinite")
  expect_error(.np_udens_q6_kernel(1,"gaussian"), "uniform or Epanechnikov")
  expect_error(.np_udens_q6_intervals(list(),c(1,Inf),.05), "maximum calibration")
})
