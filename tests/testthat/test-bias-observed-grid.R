test_that("bias-corrected panels omit unused levels without recoding the fit", {
  if (!spawn_mpi_slaves(n=1L)) skip("MPI pool unavailable in this test context")
  on.exit(close_mpi_slaves(force=TRUE),add=TRUE)
  old <- options(np.messages=FALSE);on.exit(options(old),add=TRUE)
  set.seed(505);n <- 150L;x <- runif(n)
  u <- factor(sample(c("a","b","c"),n,TRUE),levels=c("a","b","c","d"))
  d <- data.frame(x=x,u=u);y <- sin(2*pi*x)+.4*as.integer(u)+rnorm(n,sd=.25)
  b <- npregbw(xdat=d,ydat=y,bws=c(.1,.3),regtype="ll",bandwidth.compute=FALSE)
  fit <- npreg(b,txdat=d,tydat=y);saved <- fit
  v <- plot(fit,errors="asymptotic",center="bias-corrected",band="all",
    output="data",neval=2L,gaussian.draws=31L)
  expect_identical(levels(plot_coordinates(v, 2, "X")$u),levels(u))
  expect_setequal(as.character(plot_coordinates(v, 2, "X")$u),c("a","b","c"))
  expect_true(all(is.finite(plot_target(v, 2, band = "pointwise")$bias_corrected)))
  expect_identical(fit,saved)
})
