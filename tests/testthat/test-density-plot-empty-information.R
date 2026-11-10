test_that("density plot admission separates empty information from fixed-support zeros", {
  admit <- getFromNamespace(".np_plot_density_marginal_se", 'npRmpi')
  for (kernel in c("gaussian","epanechnikov","uniform","beta")) {
    b <- list(ncon=1L,ckertype=kernel,ckerorder=2L,ckerbound="fixed",
              ckerlb=0,ckerub=1,icon=TRUE)
    grid <- data.frame(x=c(-.01,0,.5,1,1.01))
    out <- suppressWarnings(admit(b,grid,rep(0,5),rep(0,5)))
    expect_identical(out,c(0,NA_real_,NA_real_,NA_real_,0),info=kernel)
    b$ckerbound <- "range"
    expect_true(all(is.na(suppressWarnings(admit(b,grid,rep(0,5),rep(0,5))))))
    b$ckerorder <- 4L
    expect_identical(admit(b,grid,rep(0,5),rep(0,5)),rep(0,5))
    b$ckerorder <- 2L
    expect_identical(admit(b,grid,rep(.2,5),rep(.01,5)),rep(.01,5))
  }
})

test_that("empty density windows retain their grid with undefined uncertainty", {
  if (!spawn_mpi_slaves(1L)) skip("MPI slaves unavailable"); on.exit(close_mpi_slaves(),add=TRUE)
  withr::local_options(np.messages=FALSE)
  x <- data.frame(x=c(seq(.01,.19,length.out=25),seq(.81,.99,length.out=25)))
  b <- npudensbw(dat=x,bws=.08,ckertype="uniform",bandwidth.compute=FALSE)
  f <- npudens(bws=b,tdat=x); saved <- f
  for (alpha in c(.05,.1)) for (band in c("pointwise","pmzsd","bonferroni","simultaneous","all")) {
    p <- suppressWarnings(plot(f,neval=21,errors="asymptotic",band=band,
                               alpha=alpha,output="data"))
    bands <- if(band=="all")c("pointwise","bonferroni","simultaneous") else band
    for (which in bands) {
      d <- as.data.frame(p,panel=1,band=which)
      empty <- vapply(d$eval_1,function(z) all(abs(x$x-z)>.08),logical(1))
      expect_true(any(empty));expect_equal(nrow(d),21L)
      expect_true(all(d$estimate[empty]==0))
      expect_true(all(is.na(d$std_error[empty])))
      expect_true(all(d$interval_status[empty]=="undefined"))
      if(which=="simultaneous") expect_true(all(is.na(d$lower))) else {
        expect_true(all(is.finite(d$lower[!empty])))
        expect_true(all(is.finite(d$upper[!empty])))
      }
    }
  }
  expect_identical(f,saved)
  raw <- npudens(bws=b,tdat=x,edat=data.frame(x=.5),se=TRUE)
  expect_identical(as.numeric(raw$dens),0)
  expect_identical(as.numeric(raw$derr),0)
  grDevices::pdf(NULL); on.exit(grDevices::dev.off(),add=TRUE)
  out <- suppressWarnings(plot(f,neval=21,errors="asymptotic",band="pointwise",output="plot-data"))
  expect_equal(nrow(as.data.frame(out,panel=1)),21L)
})
