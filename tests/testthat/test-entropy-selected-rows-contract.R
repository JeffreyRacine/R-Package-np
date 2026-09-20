test_that("selected entropy rows retain the full-sample statistic", {
  set.seed(331)
  x <- rnorm(41); y <- sin(x)+rnorm(41); bw <- c(.4,.5,.6,.7)
  call <- function(rows) .Call("C_np_entropy_bivariate_summation_rows",
    x, y, bw, rows, PACKAGE="npRmpi")
  all <- call(seq_along(x))
  scalar <- .Call("C_np_entropy_bivariate_summation", x,y,bw,PACKAGE="npRmpi")
  expect_equal(.5 * mean(all), scalar, tolerance=2e-12)
  for (rows in list(integer(),1L,41L,c(9L,2L,9L),rev(seq_along(x))))
    expect_identical(call(rows), all[rows])
  for (bad in list(NULL,NA_integer_,0L,42L,-1L,1.5))
    expect_error(call(bad), "entropy evaluation row")
  # Independent density-sum algebra, not another call to the same C owner.
  sums <- function(z,h) rowSums(exp(-.5*(outer(z,z,"-")/h)^2))
  joint <- rowSums(exp(-.5*((outer(x,x,"-")/bw[3])^2+
                           (outer(y,y,"-")/bw[4])^2)))
  oracle <- (1-sqrt(sums(x,bw[1])*sums(y,bw[2])*prod(bw[3:4])/
                   (length(x)*joint*prod(bw[1:2]))))^2
  expect_equal(all,oracle,tolerance=2e-12)
})
