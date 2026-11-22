test_that("uniform degree-zero NN means preserve explicit methods and derivative refusal", {
  skip_on_cran()
  opt <- options(np.messages=FALSE,np.plot.progress=FALSE,np.tree=FALSE)
  on.exit(options(opt),add=TRUE)
  set.seed(1126201);seed <- .Random.seed
  for(bound in c("none","fixed")) for(type in c("generalized_nn","adaptive_nn")) for(kind in c("lc","lp0","pooled")) {
    x <- data.frame(x=seq(-.97,.96,length.out=90L));y <- sin(3*x$x)+.1*cos(seq_len(90L))
    if(kind=="pooled"){x$z <- factor(rep(letters[1:3],length.out=90L));y <- y+.2*as.integer(x$z)}
    args <- list(xdat=x,ydat=y,bws=c(25,if(kind=="pooled").15),regtype=if(kind=="lp0")"lp"else"lc",bwtype=type,ckertype="uniform",ckerbound=bound,ckerlb=if(bound=="fixed")-1 else NULL,ckerub=if(bound=="fixed")1 else NULL,bandwidth.compute=FALSE)
    if(kind=="lp0")args$degree <- 0L
    bw <- do.call(npregbw,args);fit <- npreg(bw,txdat=x,tydat=y);held <- fit
    for(law in c("gaussian","mammen","rademacher","inid","fixed","geom")) {
      st <- .np_bias_regression_prepare(bw,x,y,adapt=TRUE,gradient.order=1L,law=law,calibration="absolute")
      expect_identical(st$settings$kernel.admission$type,type)
      expect_identical(st$settings$kernel.admission$kernel,"uniform")
      expect_identical(st$settings$kernel.admission$bound,bound)
      expect_false(isTRUE(st$settings$adaptation$active))
      expect_false(isTRUE(st$settings$mixed$active))
    }
    expect_true(.np_bias_regression_nn_zero(bw,"uniform"))
    expect_false(.np_bias_regression_nn_zero(bw))
    expect_error(.np_bias_regression_rows(st,x[1:3,,drop=FALSE],data.frame(row=1:3,target="d1")),"continuous derivatives remain unqualified")
    for(law in list(NULL,NA_character_,c("inid","geom"),"unknown"))
      expect_error(.np_bias_regression_prepare(bw,x,y,law=law,calibration="basic-local-envelope"),"degree-zero kernel NN inference requires")
    run <- function(law="wild",wild="mammen",cal=NULL,errors="bootstrap",band="all") {
      a <- list(x=fit,errors=errors,center="bias-corrected",neval=3L,
        gradients=FALSE,alpha=.2,random.seed=1126202L,output="data",band=band)
      if(errors=="bootstrap"){a$bootstrap <- law;a$B <- 63L}else a$gaussian.draws <- 63L
      if (!is.null(cal)) a$interval <- regression_test_interval(cal)
      if (wild != "mammen") a$boot.control <- np_boot_control(wild=wild)
      suppressWarnings(capture_plot_provider(do.call(plot,a)))
    }
    for(wild in c("mammen","rademacher")) {
      effective <- if(wild=="mammen")"basic-local-envelope"else"tail-envelope"
      default <- run(wild=wild);explicit <- run(wild=wild,cal=effective)
      a <- default$provider[[1L]]$inference;e <- explicit$provider[[1L]]$inference
      expect_identical(a$calibration.method,effective);expect_identical(a$calibration.requested,"default")
      a$calibration.requested <- e$calibration.requested
      expect_identical(a,e);expect_identical(default$result$panels,explicit$result$panels)
      value <- run(wild=wild,cal="symmetric-sd")
      expect_s3_class(value$result,"np_plot_data")
      expect_identical(value$provider[[1L]]$inference$calibration.method,"absolute")
    }
    expect_s3_class(run(errors="asymptotic")$result,"np_plot_data")
    expect_s3_class(run(band="pointwise",cal="normal-sd")$result,"np_plot_data")
    expect_s3_class(run(law="inid",cal="absolute")$result,"np_plot_data")
    expect_s3_class(run(law="inid",cal="basic-envelope")$result,"np_plot_data")
    expect_error(run(cal="basic-envelope"),"basic-envelope requires inid, fixed or geom")
    higher.args <- args;higher.args$regtype <- "lp";higher.args$degree <- 1L
    higher <- do.call(npregbw,higher.args)
    expect_identical(.np_bias_regression_kernel_class(higher,"absolute"),"unsupported")
    expect_error(.np_bias_regression_prepare(higher,x,y,law="mammen",calibration="basic-local-envelope"),"requires a supported second-order kernel")
    expect_identical(.Random.seed,seed);expect_identical(fit,held)
  }
})

test_that("uniform NN main rows preserve physical support ties", {
  skip_on_cran()
  opt <- options(np.messages=FALSE,np.plot.progress=FALSE,np.tree=FALSE)
  on.exit(options(opt),add=TRUE)
  set.seed(4202)
  x <- data.frame(x=c(18,67,sample(18:67,148L,replace=TRUE)))
  y <- log(x$x)+rnorm(150L,sd=.1)
  seed <- .Random.seed
  for(type in c("generalized_nn","adaptive_nn")) for(bound in c("none","fixed")) for(reg in c("lc","lp")) {
    bw <- suppressWarnings(npregbw(xdat=x,ydat=y,bws=30,bwtype=type,
      ckertype="uniform",ckerbound=bound,ckerlb=if(bound=="fixed")18 else NULL,
      ckerub=if(bound=="fixed")67 else NULL,regtype=reg,degree=0L,bandwidth.compute=FALSE))
    fit <- npreg(bw,txdat=x,tydat=y); held <- fit
    value <- suppressWarnings(capture_plot_provider(plot(fit,center="bias-corrected",
      errors="bootstrap",bootstrap="inid",neval=50L,B=31L,
      random.seed=1120501L,output="data")))
    inf <- value$provider[[1L]]$inference
    grid <- inf$grid
    # Independent physical-coordinate LC weights, including donor-specific
    # bounded ANN mass. The strict support convention is unchanged.
    cdf <- function(z).5+.5*pmax(-1,pmin(1,z))
    radius <- if(type=="adaptive_nn") vapply(seq_len(nrow(x)),function(i)
      sort(abs(x$x[-i]-x$x[i]))[30L],0) else NULL
    expected <- vapply(grid$x,function(g) {
      h <- if(type=="generalized_nn")sort(abs(x$x-g))[30L] else radius
      at <- if(type=="generalized_nn")g else x$x
      mass <- if(bound=="fixed")cdf((67-at)/h)-cdf((18-at)/h) else 1
      w <- .5*(abs((x$x-g)/h)<1)/h/mass
      sum(w*y)/sum(w)
    },0)
    expect_equal(inf$estimate,expected,tolerance=1e-12)
    expect_equal(inf$estimate,npreg(bw,txdat=x,tydat=y,exdat=grid)$mean,tolerance=1e-12)
    expect_identical(fit,held);expect_identical(.Random.seed,seed)
  }
})
