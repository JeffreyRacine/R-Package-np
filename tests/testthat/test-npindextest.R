test_that("index test is exported with registered methods and validates arguments", {
  ns <- asNamespace("npRmpi")
  expect_identical(getExportedValue("npRmpi","npindextest"), ns$npindextest)
  for (generic in c("print","summary"))
    expect_identical(getS3method(generic,"indextest",envir=ns),
                     get(paste0(generic,".indextest"),ns,inherits=FALSE))
  expect_error(npindextest(NULL,B=2L),"at least 9")
  expect_error(npindextest(NULL),"object must")
  expect_error(npindextest(NULL,bw.method="cv"),"arg")
})

test_that("the native index score agrees with the ordered-pair statistic", {
  x <- cbind(c(-1, 0, 1, .5, -.5), c(1, .2, -.5, 1.2, 1.4))
  y <- c(.5, 1, 1.5, -1, 2)
  z <- x[,1]
  g <- 0.4
  h <- c(.3,.5)
  dz <- outer(z,z,"-")
  dist2 <- as.double(dist(x))^2
  u <- rowSums(outer(y,y,"-")*exp(-.5*(dz/g)^2))
  base <- tcrossprod(u)*exp(-.5*outer(x[,2],x[,2],"-")^2)
  diag(base) <- 0
  ref <- vapply(h,function(width) {
    a <- base*exp(-.5*(dz/width)^2)
    sqrt(5/4)*sum(a)/sqrt(2*sum(a*a))
  },numeric(1))
  expect_equal(.Call(C_np_indextest_score,y,z,dist2,g,h),ref,tolerance=1e-12)
})
