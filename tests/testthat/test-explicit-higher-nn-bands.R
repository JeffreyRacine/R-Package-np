test_that("explicit higher Epanechnikov NN bands retain the fitted hat and the requested target", {
  skip_on_cran()
  opt <- options(np.messages=FALSE, np.plot.progress=FALSE)
  on.exit(options(opt),add=TRUE)
  x <- data.frame(x=seq(-.9,.98,length.out=100))
  y <- sin(3*x$x)+.1*cos(seq_len(nrow(x)))
  b <- npregbw(xdat=x,ydat=y,bws=26,regtype="lp",degree=3L,basis="glp",
    bernstein.basis=TRUE,bwtype="generalized_nn",ckertype="epanechnikov",bandwidth.compute=FALSE)
  fit <- npreg(b,txdat=x,tydat=y)
  a <- list(x=fit,errors="bootstrap",center="bias-corrected",neval=3L,B=99L,
    random.seed=146L,alpha=.1,output="data",gradients=TRUE,band="simultaneous")
  expect_error(do.call(plot,a))
  for (construction in c("symmetric-sd","basic-local-envelope")) {
    z <- suppressWarnings(do.call(plot,c(a,list(interval=construction))))
    m <- z$inference[[1L]]
    expect_identical(m$resampling_law,"mammen")
    expect_identical(m$components$simultaneous$interval,construction)
  }
  expect_true(.np_bias_explicit_nn_target(utils::modifyList(b,list(type="adaptive_nn")),TRUE))
  expect_false(.np_bias_explicit_nn_target(utils::modifyList(b,list(type="adaptive_nn",degree=4L)),TRUE))
  expect_false(.np_bias_explicit_nn_target(utils::modifyList(b,list(bernstein.basis=FALSE)),FALSE))
  st <- .np_bias_regression_prepare(b,x,y,adapt=TRUE,gradient.order=1L,law="mammen",
    calibration="absolute",explicit.nn=TRUE)
  g <- x[round(seq(1, nrow(x), length.out=25L)),,drop=FALSE]
  rows <- .np_bias_regression_rows(st,g,data.frame(row=seq_len(nrow(g)),target="d1"))
  H <- npreghat(b,txdat=x,exdat=g,s=1L,output="matrix")
  expect_equal(as.double(rows$original),as.double(H),tolerance=1e-10)
  expect_error(.np_bias_main_admit(st,matrix(1),1e-12,list(degree.engine=3L)),"zero-ridge")
  expect_error(.np_bias_regression_rows(st,data.frame(x=1),data.frame(row=1L,target="d1")),"observed training range")
})
test_that("positive-degree beta NN slope choices are explicit", {
  skip_on_cran()
  options0 <- options(np.messages=FALSE,np.plot.progress=FALSE)
  on.exit(options(options0),add=TRUE)
  x <- data.frame(x=seq(-.9,.95,length.out=100L))
  y <- exp(x$x)+.05*cos(seq_len(nrow(x)))
  for (type in c("generalized_nn","adaptive_nn")) for (degree in if(type=="generalized_nn") 1:2 else 2L) {
    b <- npregbw(xdat=x,ydat=y,bws=29,regtype="lp",degree=degree,basis="glp",bernstein.basis=TRUE,
                 bwtype=type,ckertype="beta",ckerbound="fixed",ckerlb=-1,ckerub=1,bandwidth.compute=FALSE)
    a <- list(x=npreg(b,txdat=x,tydat=y),gradients=TRUE,neval=3L,center="bias-corrected",errors="bootstrap",
              output="data",band="pointwise",B=99L,alpha=.1,random.seed=31L)
    expect_error(do.call(plot,a))
    for (interval in c("symmetric-sd","basic-local-envelope")) {
      z <- suppressWarnings(do.call(plot,c(a,list(interval=interval))))
      expect_identical(z$inference[[1L]]$components$pointwise$interval,interval)
    }
    if (type == "generalized_nn" && degree == 2L) {
      a$band <- "simultaneous"
      absolute <- suppressWarnings(do.call(plot,c(a,list(interval="symmetric-sd"))))
      expect_identical(absolute$inference[[1L]]$components$simultaneous$interval,"symmetric-sd")
      for (scope in c("bonferroni", "all")) {
        companion <- a
        companion$band <- scope
        expect_error(do.call(plot,c(companion,list(interval="symmetric-sd"))), "mean target")
      }
      local <- suppressWarnings(do.call(plot,c(a,list(interval="basic-local-envelope"))))
      expect_identical(local$inference[[1L]]$components$simultaneous$interval,"basic-local-envelope")
    }
  }
})
test_that("experimental linear uniform ANN mean bands keep continuous slopes outside admission", {
  skip_on_cran()
  old <- options(np.messages=FALSE,np.plot.progress=FALSE)
  on.exit(options(old),add=TRUE)
  x <- data.frame(x=seq(.02,.98,length.out=90))
  y <- x$x^2+.05*cos(seq_len(nrow(x)))
  b <- npregbw(xdat=x,ydat=y,bws=19,regtype="lp",degree=1L,basis="glp",bernstein.basis=TRUE,
               bwtype="adaptive_nn",ckertype="uniform",bandwidth.compute=FALSE)
  args <- list(x=npreg(b,txdat=x,tydat=y),errors="bootstrap",center="bias-corrected",band="simultaneous",
               neval=3L,alpha=.1,B=99L,random.seed=23L,output="data")
  expect_error(do.call(plot,args))
  args$interval <- "basic-local-envelope"
  r <- suppressWarnings(do.call(plot,args))
  expect_identical(r$inference[[1L]]$resampling_law,"mammen")
  expect_identical(r$inference[[1L]]$components$simultaneous$interval,args$interval)
  expect_error(do.call(plot,c(args,list(gradients=TRUE))))
  expect_false(.np_bias_explicit_nn_target(utils::modifyList(b,list(ckerbound="fixed")),FALSE))
})

test_that("explicit Epanechnikov ANN first derivatives retain their row guards", {
  skip_on_cran()
  old <- options(np.messages=FALSE, np.plot.progress=FALSE)
  on.exit(options(old), add=TRUE)
  v <- readRDS(test_path("fixtures", "r112-epanechnikov-ann.rds"))
  x <- data.frame(x=v$x)
  b <- npregbw(xdat=x, ydat=v$y, bws=53, bandwidth.compute=FALSE,
    regtype="lp", degree=3L, basis="glp", bernstein.basis=TRUE,
    bwtype="adaptive_nn", ckertype="epanechnikov")
  fit <- npreg(b, txdat=x, tydat=v$y)
  a <- list(x=fit, errors="bootstrap", center="bias-corrected", band="simultaneous",
    gradients=TRUE, output="data", B=99L, neval=3L, random.seed=91L)
  expect_error(do.call(plot, a))
  got <- suppressWarnings(do.call(plot, c(a, list(interval="symmetric-sd"))))
  expect_identical(got$inference[[1L]]$components$simultaneous$interval,"symmetric-sd")
  st <- .np_bias_regression_prepare(b, x, v$y, adapt=TRUE, gradient.order=1L,
    law="mammen", calibration="absolute", explicit.nn=TRUE)
  outside <- data.frame(x=1.1496080229039194)
  H <- npreghat(b, txdat=x, exdat=outside, s=1L)
  rho <- attr(H, "ridge.used", exact=TRUE)
  expect_gt(rho, 0)
  expect_error(.np_bias_regression_rows(st, outside, data.frame(row=1L,target="d1")),
    "observed training range")
  spec <- st
  spec$x <- st$main.physical$x
  spec$bws <- st$main.physical$bws
  expect_error(.np_bias_main_rows(spec, outside, 1L), "identifiable full-rank")
  expect_error(.np_bias_main_admit(spec, H, rho, list(degree.engine=3L)), "zero-ridge")
})

test_that("cubic beta ANN bands refuse the whole grid when an endpoint selects ridge", {
  skip_on_cran()
  opt <- options(np.messages=FALSE, np.plot.progress=FALSE)
  on.exit(options(opt), add=TRUE)
  v <- readRDS(test_path("fixtures", "r112-cubic-beta-ann.rds"))
  x <- data.frame(x=v$x)
  bw <- npregbw(xdat=x, ydat=v$y, bws=24L, bandwidth.compute=FALSE,
    regtype="lp", degree=3L, basis="glp", bernstein.basis=TRUE,
    bwtype="adaptive_nn", ckertype="beta", ckerbound="fixed", ckerlb=-1, ckerub=1)
  fit <- npreg(bw, txdat=x, tydat=v$y)
  a <- list(x=fit, errors="bootstrap", center="bias-corrected", B=59L, neval=3L,
            random.seed=97L, output="data", interval="basic-local-envelope")
  for (target in c(FALSE, TRUE))
    expect_error(do.call(plot, c(a, list(gradients=target, band="simultaneous"))), "zero-ridge")
  expect_error(do.call(plot, c(a, list(band="all"))))
  expect_error(do.call(plot, c(a, list(band="pointwise"))))
})

test_that("known-support quadratic Gaussian GNN uses explicit Mammen local bands", {
  skip_on_cran()
  opt <- options(np.messages=FALSE, np.plot.progress=FALSE)
  on.exit(options(opt), add=TRUE)
  x <- data.frame(x=seq(-.95,.95,length.out=80L))
  y <- sin(3*x$x)+.1*cos(seq_len(nrow(x)))
  bw <- npregbw(xdat=x, ydat=y, bws=18, bandwidth.compute=FALSE, regtype="lp",
    degree=2L, basis="glp", bernstein.basis=TRUE, bwtype="generalized_nn",
    ckertype="gaussian", ckerbound="fixed", ckerlb=-1, ckerub=1)
  a <- list(x=npreg(bw,txdat=x,tydat=y), center="bias-corrected", errors="bootstrap",
            gradients=TRUE, neval=4L, B=79L, output="data", band="all", random.seed=32L)
  expect_error(do.call(plot,a))
  expect_error(do.call(plot,c(a,list(interval="symmetric-sd"))))
  z <- suppressWarnings(do.call(plot,c(a,list(interval="basic-local-envelope"))))
  expect_identical(z$inference[[1L]]$resampling_law,"mammen")
  expect_identical(z$inference[[1L]]$components$bonferroni$interval,"basic-local-envelope")
  pairs <- suppressWarnings(do.call(plot,c(a,list(interval="basic-local-envelope",bootstrap="inid"))))
  expect_identical(pairs$inference[[1L]]$resampling_law,"inid")
  mean.args <- a
  mean.args$gradients <- FALSE
  absolute <- suppressWarnings(do.call(plot,c(mean.args,list(interval="symmetric-sd"))))
  expect_identical(absolute$inference[[1L]]$components$simultaneous$interval,"symmetric-sd")
  expect_error(do.call(plot,mean.args))
  adaptive <- utils::modifyList(bw,list(type="adaptive_nn"))
  expect_true(.np_bias_explicit_nn_target(adaptive,TRUE))
  expect_false(.np_bias_gnn_local_pairs(adaptive,TRUE,"inid","basic-local-envelope"))
  for (fields in list(list(ncon=2L),list(ckerub=Inf),
                      list(degree=3L),list(bernstein.basis=FALSE)))
    expect_false(.np_bias_explicit_nn_target(utils::modifyList(bw,fields),TRUE))
})

test_that("explicit local Epanechnikov GNN pairs retain the cubic slope boundary", {
  skip_on_cran()
  opt <- options(np.messages=FALSE, np.plot.progress=FALSE)
  on.exit(options(opt),add=TRUE)
  x <- data.frame(x=seq(-.89,.97,length.out=80L))
  y <- sin(x$x)+.03*cos(seq_len(nrow(x)))
  b <- npregbw(xdat=x,ydat=y,bws=19,regtype="lp",degree=2L,basis="glp",
    bernstein.basis=TRUE,bwtype="generalized_nn",ckertype="epanechnikov",bandwidth.compute=FALSE)
  fit <- npreg(b,txdat=x,tydat=y)
  a <- list(x=fit,errors="bootstrap",center="bias-corrected",gradients=TRUE,
            bootstrap="fixed",boot.control=np_boot_control(blocklen=5L),output="data",
            neval=4L,B=69L,band="all",random.seed=231L)
  expect_error(do.call(plot,a))
  obj <- suppressWarnings(do.call(plot,c(a,list(interval="basic-local-envelope"))))
  expect_identical(obj$inference[[1L]]$resampling_law,"fixed")
  expect_identical(obj$inference[[1L]]$components$bonferroni$interval,"basic-local-envelope")
  expect_false(.np_bias_gnn_local_pairs(utils::modifyList(b,list(degree=3L)),TRUE,"inid","basic-local-envelope"))
  expect_false(.np_bias_gnn_local_pairs(b,TRUE,"rademacher","basic-local-envelope"))
})
