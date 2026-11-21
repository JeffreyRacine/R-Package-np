test_that("Epanechnikov first-order NN preserves defaults and partial derivatives", {
  opt <- options(np.messages=FALSE,np.plot.progress=FALSE,np.tree=FALSE)
  on.exit(options(opt),add=TRUE)
  set.seed(1127601);seed <- .Random.seed
  for(type in c("generalized_nn","adaptive_nn"))for(bound in c("none","fixed"))for(kind in c("q1","pooled","partial")) {
    n <- 96L
    x <- data.frame(x=seq(-.97,.96,length.out=n))
    y <- sin(3*x$x)+.1*cos(seq_len(n))
    if(kind=="pooled")x$z <- factor(rep(letters[1:3],length.out=n))
    if(kind=="partial")x$x2 <- .9*sin(seq_len(n)*sqrt(2))
    nc <- which(vapply(x,is.numeric,logical(1)))
    degree <- if(kind=="partial")c(0L,1L)else 1L
    bw <- npregbw(xdat=x,ydat=y,bws=if(kind=="partial")c(45,45)else c(25,if(kind=="pooled").15),
      regtype="lp",degree=degree,basis="glp",bernstein.basis=TRUE,bwtype=type,
      ckertype="epanechnikov",ckerbound=bound,ckerlb=if(bound=="fixed")rep(-1,length(nc))else NULL,ckerub=if(bound=="fixed")rep(1,length(nc))else NULL,
      bandwidth.compute=FALSE)
    fit <- suppressWarnings(npreg(bw,txdat=x,tydat=y,gradients=TRUE));held <- fit
    run <- function(law="mammen",cal=NULL,errors="bootstrap",gradients=FALSE) {
      a <- list(x=fit,errors=errors,center="bias-corrected",neval=3L,
        gradients=gradients,alpha=.2,random.seed=1127602L,output="data",band="all")
      if(errors=="bootstrap"){a$bootstrap <- if(law%in%c("mammen","rademacher"))"wild"else law;a$B <- 63L}else a$gaussian.draws <- 63L
      if(!is.null(cal)||law=="rademacher")a$boot.control <- np_boot_control(calibration=if(is.null(cal))"default"else cal,wild=if(law=="rademacher")"rademacher"else"mammen")
      suppressWarnings(capture_plot_provider(do.call(plot,a)))
    }
    for(law in c("mammen","inid","rademacher")) {
      effective <- if(law=="rademacher")"tail-envelope"else"basic-local-envelope"
      a <- run(law);b <- run(law,effective)
      expect_identical(a$result$panels,b$result$panels)
      for(j in seq_along(a$provider)) {
        ai <- a$provider[[j]]$inference;bi <- b$provider[[j]]$inference
        expect_identical(ai$calibration.method,effective)
        expect_identical(ai$calibration.requested,"default")
        ai$calibration.requested <- bi$calibration.requested
        expect_identical(ai,bi)
      }
    }
    expect_s3_class(run(errors="asymptotic")$result,"np_plot_data")
    grad <- run(gradients=TRUE)
    expect_s3_class(grad$result,"np_plot_data")
    if(kind=="partial") {
      expect_true(all(is.na(grad$provider[[1L]]$grad)))
      expect_true(any(is.finite(grad$provider[[2L]]$grad)))
      expect_error(npreghat(bw,txdat=x,exdat=x[1:3,,drop=FALSE],s=c(1L,0L),output="matrix"),
        "exceeds local polynomial degree")
    }
    expect_s3_class(run(law="inid",cal="basic-envelope")$result,"np_plot_data")
    expect_s3_class(run(cal="symmetric-sd")$result,"np_plot_data")
    expect_error(run(cal="basic-envelope"),"basic-envelope requires")
    for(law in list(NULL,NA_character_,c("inid","geom"),"unknown"))
      expect_error(.np_bias_regression_prepare(bw,x,y,law=law,calibration="basic-local-envelope"),"first-order Epanechnikov NN inference requires")
    outside <- bw;outside$ckertype <- "uniform"
    expect_identical(.np_bias_regression_kernel_class(outside,"absolute"),"unsupported")
    outside$ckertype <- "beta"
    expect_identical(.np_bias_regression_kernel_class(outside,"absolute"),"unsupported")
    outside <- bw;outside$degree <- rep(2L,length(nc))
    expect_identical(.np_bias_regression_kernel_class(outside,"absolute"),"unsupported")
    expect_identical(fit,held);expect_identical(.Random.seed,seed)
  }
})

test_that("Epanechnikov local linear agrees with its raw GLP model", {
  opt <- options(np.messages=FALSE,np.plot.progress=FALSE,np.tree=FALSE)
  on.exit(options(opt),add=TRUE)
  x0 <- data.frame(x=seq(-.97,.96,length.out=96L));y <- sin(3*x0$x)+.1*cos(seq_len(96L))
  set.seed(1128201);seed <- .Random.seed
  for(type in c("generalized_nn","adaptive_nn"))for(mixed in c(FALSE,TRUE)) {
    x <- x0;if(mixed)x$z <- factor(rep(letters[1:3],length.out=96L))
    args <- list(xdat=x,ydat=y,bws=c(25,if(mixed).15),bwtype=type,ckertype="epanechnikov",
      ckerbound="fixed",ckerlb=-1,ckerub=1,bandwidth.compute=FALSE)
    ll <- do.call(npregbw,c(args,list(regtype="ll")))
    lp <- do.call(npregbw,c(args,list(regtype="lp",degree=1L,basis="glp",bernstein.basis=FALSE)))
    fit.ll <- npreg(ll,txdat=x,tydat=y);fit.lp <- npreg(lp,txdat=x,tydat=y)
    for(gradient in c(FALSE,TRUE)) {
      a <- suppressWarnings(plot(fit.ll,gradients=gradient,errors="bootstrap",center="bias-corrected",
        bootstrap="wild",B=63L,random.seed=1128202L,neval=3L,alpha=.2,band="all",output="data"))
      b <- suppressWarnings(plot(fit.lp,gradients=gradient,errors="bootstrap",center="bias-corrected",
        bootstrap="wild",B=63L,random.seed=1128202L,neval=3L,alpha=.2,band="all",output="data"))
      expect_identical(a$panels,b$panels)
    }
    expect_identical(.Random.seed,seed)
  }
})
