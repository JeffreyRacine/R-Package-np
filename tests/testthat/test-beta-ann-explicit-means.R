test_that("linear bounded beta ANN mean intervals require an explicit assessed construction", {
  skip_on_cran()
  if (!spawn_mpi_slaves(n=1L)) skip("MPI workers unavailable")
  on.exit(close_mpi_slaves(force=TRUE),add=TRUE)
  opt <- options(np.messages=FALSE,np.plot.progress=FALSE)
  on.exit(options(opt),add=TRUE)
  set.seed(73011L)
  x <- data.frame(x=seq(-.93,.98,length.out=100L))
  y <- sin(3*x$x)+.1*cos(seq_len(nrow(x)))
  b <- npregbw(xdat=x,ydat=y,bws=34,regtype="lp",degree=1L,basis="glp",
    bernstein.basis=TRUE,bwtype="adaptive_nn",ckertype="beta",
    ckerbound="fixed",ckerlb=-1,ckerub=1,bandwidth.compute=FALSE)
  fit <- npreg(b,txdat=x,tydat=y)
  a <- list(x=fit, errors="bootstrap",center="bias-corrected",neval=3L,
            B=99L,alpha=.1,random.seed=146L,output="data",band="simultaneous")
  held <- .Random.seed
  for (scope in c("simultaneous","pointwise")) {
    a$band <- scope
    for (construction in c("symmetric-sd","basic-local-envelope")) {
      z <- suppressWarnings(do.call(plot,c(a,list(interval=construction))))
      meta <- z$inference[[1]]
      expect_identical(meta$resampling_law,"mammen")
      expect_identical(meta$components[[scope]]$interval,construction)
      endpoints <- as.data.frame(z,band=scope)
      expect_true(all(is.finite(endpoints$lower) & endpoints$lower <= endpoints$upper))
    }
  }
  expect_error(do.call(plot,a),"supported second-order")
  expect_error(do.call(plot,c(a,list(boot.control=np_boot_control(wild="mammen")))),"supported second-order")
  a$interval <- "symmetric-sd"
  for (s in c("bonferroni","all")) {
    result <- suppressWarnings(do.call(plot,utils::modifyList(a,list(band=s))))
    expect_identical(result$inference[[1L]]$components$bonferroni$interval,"symmetric-sd")
  }
  for (extra in list(list(bootstrap="inid"),
                     list(boot.control=np_boot_control(wild="rademacher"))))
    expect_error(do.call(plot,utils::modifyList(a,extra)),"supported second-order")
  expect_error(do.call(plot,utils::modifyList(a,list(bias.control=np_bias_control(iterations=1L)))),
               "auxiliary bias settings")
  slope <- utils::modifyList(a,list(gradients=TRUE))
  for (scope in c("simultaneous","pointwise")) {
    slope$band <- scope
    for (interval in c("symmetric-sd","basic-local-envelope")) {
      slope$interval <- interval
      obj <- suppressWarnings(do.call(plot,slope))
      expect_identical(obj$inference[[1L]]$components[[scope]]$interval,interval)
      limits <- as.data.frame(obj,band=scope)
      expect_true(all(is.finite(limits$lower) & limits$lower <= limits$upper))
    }
  }
  slope$interval <- NULL
  expect_error(do.call(plot,slope),"supported second-order")
  slope$interval <- "symmetric-sd"
  result <- suppressWarnings(do.call(plot,utils::modifyList(slope,list(band="all"))))
  expect_identical(result$inference[[1L]]$components$bonferroni$interval,"symmetric-sd")
  expect_true(.np_bias_beta_nn_mean(b,"adaptive_nn",1L))
  for (change in list(list(bernstein.basis=FALSE),list(degree=2L),list(ncon=2L),
                      list(nuno=1L),list(ckerub=Inf),list(ckerbound="none")))
    expect_false(.np_bias_beta_nn_mean(utils::modifyList(b,change),"adaptive_nn",1L))
  st <- .np_bias_regression_prepare(b,x,y,law="mammen",calibration="absolute",explicit.ann.mean=TRUE)
  expect_error(.np_bias_regression_rows(st,data.frame(x=-.99),data.frame(target="mean")),
               "observed training range")
  expect_error(.np_bias_main_admit(st,matrix(1),1e-12,list(degree.engine=1L)),"zero-ridge")
  expect_identical(.Random.seed,held)
})
