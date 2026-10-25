test_that("bias-corrected panels omit unused levels without recoding the fit", {
  old <- options(np.messages=FALSE);on.exit(options(old),add=TRUE)
  set.seed(505);n <- 150L;x <- runif(n)
  u <- factor(sample(c("a","b","c"),n,TRUE),levels=c("a","b","c","d"))
  d <- data.frame(x=x,u=u);y <- sin(2*pi*x)+.4*as.integer(u)+rnorm(n,sd=.25)
  b <- npregbw(xdat=d,ydat=y,bws=c(.1,.3),regtype="ll",bandwidth.compute=FALSE)
  fit <- npreg(b,txdat=d,tydat=y);saved <- fit
  v <- plot(fit,errors="asymptotic",center="bias-corrected",band="all",
    output="data",neval=2L,gaussian.draws=31L)
  expect_identical(levels(v[[2]]$inference$grid$u),levels(u))
  expect_setequal(as.character(v[[2]]$inference$grid$u),c("a","b","c"))
  expect_true(all(is.finite(v[[2]]$inference$center)))
  expect_identical(fit,saved)
})
