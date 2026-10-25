test_that("categorical correction reduces to cell means and residual wild covariance", {
  if (!spawn_mpi_slaves(n=1L)) skip("MPI pool unavailable in this test context")
  on.exit(close_mpi_slaves(force=TRUE),add=TRUE)
  if (!spawn_mpi_slaves(n=1L)) skip("MPI pool unavailable in this test context")
  on.exit(close_mpi_slaves(force=TRUE),add=TRUE)
  old <- options(np.messages=FALSE);on.exit(options(old),add=TRUE)
  set.seed(61006);x <- data.frame(u=factor(rep(1:3,each=20)));y <- as.integer(x$u)+rnorm(60)
  b <- npregbw(xdat=x,ydat=y,bws=.2,regtype="lc",bandwidth.compute=FALSE)
  a <- npRmpi:::.np_bias_regression_prepare(b,x,y);g <- x[c(1,21,41),,drop=FALSE]
  v <- npRmpi:::.np_bias_regression_rows(a,g,data.frame(row=1:3,target="mean"))
  expect_equal(v$center,as.numeric(tapply(y,x$u,mean)),tolerance=1e-13)
  expect_equal(v$se^2,as.numeric(tapply(y,x$u,var))/20,tolerance=1e-13)
  expect_equal(tcrossprod(sweep(v$W,2L,a$residual,"*")),tcrossprod(v$factor),tolerance=1e-13)
  for(law in c("gaussian","mammen")) {
    z <- npRmpi:::.np_bias_inference(a,g[2,,drop=FALSE],"c1",.05,63L,law)
    expect_identical(z$bands$pointwise,z$bands$simultaneous)
    expect_identical(z$bands$pointwise,z$bands$bonferroni)
    expect_true(all(is.finite(z$bands$simultaneous)))
  }
  expect_error(npRmpi:::.np_bias_regression_prepare(b,x,y,np_bias_control(p0.scale=4)),"continuous auxiliary geometry")
  expect_error(npRmpi:::.np_bias_regression_prepare(b,x,y,np_bias_control(adaptation="local-to-global")),"requires a continuous predictor")
  fixed.scale <- npRmpi:::.np_bias_inference(a,g[2,,drop=FALSE],"mean",.05,31L,"gaussian",studentize=FALSE)
  expect_identical(fixed.scale$covariance$studentization,"within-profile residual SD")
  expect_identical(a$settings$iterations,16L)
})


test_that("categorical-only public panels preserve declared levels and references", {
  if (!spawn_mpi_slaves(n=1L)) skip("MPI pool unavailable in this test context")
  on.exit(close_mpi_slaves(force=TRUE),add=TRUE)
  old <- options(np.messages=FALSE);on.exit(options(old),add=TRUE)
  set.seed(61016)
  x <- data.frame(u=factor(rep(c("a","b","c"),each=20),
    levels=c("a","b","c","unused")))
  y <- as.integer(x$u)+rnorm(nrow(x))
  b <- npregbw(xdat=x,ydat=y,bws=.2,regtype="lc",bandwidth.compute=FALSE)
  fit <- npreg(b,txdat=x,tydat=y);saved <- fit;rng <- .Random.seed
  for (gradients in c(FALSE,TRUE)) {
    z <- plot(fit,errors="asymptotic",center="bias-corrected",band="all",
      gradients=gradients,output="data",gaussian.draws=31L)[[1L]]$inference
    expect_identical(levels(z$grid$u),levels(x$u))
    expect_identical(as.character(z$grid$u),c("a","b","c"))
    expect_equal(z$center,if(gradients)as.numeric(tapply(y,x$u,mean)[1:3]-mean(y[x$u=="a"]))
      else as.numeric(tapply(y,x$u,mean)[1:3]),tolerance=1e-13)
    expect_identical(z$settings$iterations,16L)
  }
  expect_identical(fit,saved)
  expect_identical(rng,.Random.seed)
  missing <- x;missing$u <- factor(as.character(x$u),levels=c("absent","a","b","c"))
  b <- npregbw(xdat=missing,ydat=y,bws=.2,regtype="lc",bandwidth.compute=FALSE)
  expect_error(plot(b,errors="asymptotic",center="bias-corrected",band="all",
    gradients=TRUE,output="data",gaussian.draws=31L),"observed endpoint profiles")
})
