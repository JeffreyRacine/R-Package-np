test_that("P1 higher derivatives reproduce moments in concentrated weights", {
  u <- seq(-1,1,length.out=60L)
  Q <- outer(u,0:5,`^`)
  at <- .9; h <- .15
  weight <- exp(-(u-at)^2/(2*h^2))
  ratios <- np:::.np_bias_gaussian_derivatives((u-at)/h,h,3L)
  values <- t(vapply(0:3,function(j) vapply(0:5,function(k)
    if(k<j)0 else factorial(k)/factorial(k-j)*at^(k-j),numeric(1)),numeric(6)))
  row <- np:::.np_bias_pilot_higher(Q,weight,values,ratios)
  expect_equal(drop(row %*% Q),values[4L,],tolerance=1e-9)
  expect_error(np:::.np_bias_pilot_higher(cbind(Q,Q[,2L]),weight,
    cbind(values,values[,2L]),ratios),"numerically unidentified",fixed=TRUE)
  # A fixed penalty enters the zeroth Gram matrix only. Its derivative is zero.
  P <- Q[,1:3]; V <- values[,1:3]
  ordinary <- np:::.np_bias_gaussian_derivatives((u-at)/.8,.8,3L)
  a <- np:::.np_bias_pilot_higher(P,rep(1,length(u)),V,ordinary,
    penalty=rep(.1,3))
  b <- np:::.np_bias_weighted_derivative(P,rep(1,length(u)),V,ordinary,
    penalty=rep(.1,3))
  expect_equal(a,b,tolerance=1e-10)
})
