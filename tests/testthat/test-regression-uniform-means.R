test_that("uniform iterated mean admission keeps support and law boundaries", {
  b<-list(type="fixed",ckertype="uniform",ckerorder=2L,ckerbound="none",ncon=1L,regtype="lc")
  for(bound in c("none","fixed"))for(cal in c("absolute","separate-tails","tail-envelope","basic-envelope","basic-local-envelope")) {
    b$ckerbound<-bound;expect_identical(.np_bias_regression_kernel_class(b,cal),"local-envelope")
  }
  for(type in c("generalized_nn","adaptive_nn")){z<-b;z$type<-type;expect_identical(.np_bias_regression_kernel_class(z,"basic-local-envelope"),"local-envelope")}
  for(bound in c("range","lbound","ubound")){z<-b;z$ckerbound<-bound;expect_identical(.np_bias_regression_kernel_class(z,"absolute"),"unsupported")}
  for (degree in list(1L, 2L, 3L, c(0L,1L))) {
    z <- b
    z$regtype <- "lp"
    z$basis <- "glp"
    z$degree <- degree
    for (cal in c("absolute","separate-tails","tail-envelope","basic-envelope","basic-local-envelope"))
      expect_identical(.np_bias_regression_kernel_class(z,cal), "unsupported")
  }
  z<-b;z$ckerorder<-4L;expect_identical(.np_bias_regression_kernel_class(z,"absolute"),"unsupported")
  for(law in c("wild","inid","fixed","geom"))expect_identical(.np_bias_regression_calibration(b,"default","bootstrap","all",law,"mammen"),"basic-local-envelope")
  expect_identical(.np_bias_regression_calibration(b,"default","bootstrap","all","wild","rademacher"),"tail-envelope")
  expect_identical(.np_bias_regression_calibration(b,"default","bootstrap","pmzsd","wild","rademacher"),"absolute")
  expect_identical(.np_bias_regression_calibration(b,"default","asymptotic","all",NULL,NULL),"absolute")
  for(cal in c("absolute","separate-tails","tail-envelope","basic-envelope","basic-local-envelope"))expect_identical(.np_bias_regression_calibration(b,cal,"bootstrap","all","wild","mammen"),cal)
})

test_that("uniform main mean rows match an independent weighted SVD", {
  skip_on_cran()
  if(!spawn_mpi_slaves(n=1L))skip("MPI workers unavailable")
  on.exit(close_mpi_slaves(force=TRUE),add=TRUE)
  old<-options(np.messages=FALSE,np.plot.progress=FALSE);on.exit(options(old),add=TRUE)
  degree <- 0L
  for(pooled in c(FALSE,TRUE))for(bound in c("none","fixed")) {
    x<-data.frame(x=if(pooled)rep(seq(-.9,.9,length.out=30L),3L)else seq(-.9,.9,length.out=90L))
    if(pooled)x$group<-factor(rep(letters[1:3],each=30L))
    y<-sin(3*x$x)+.15*cos(seq_len(90L))+if(pooled).2*as.integer(x$group)else 0
    args<-list(xdat=x,ydat=y,bws=c(1.1,if(pooled).2),regtype="lp",degree=degree,ckertype="uniform",bandwidth.compute=FALSE,ckerbound=bound)
    if(bound=="fixed")args<-c(args,list(ckerlb=-.9,ckerub=.9))
    bw<-do.call(npregbw,args);held<-bw;grid<-x[c(3L,33L,63L,17L,47L,77L),,drop=FALSE]
    st<-.np_bias_regression_prepare(bw,x,y,adapt=TRUE,gradient.order=1L,law="mammen",calibration="basic-local-envelope")
    expect_identical(st$settings$main$kernel,"uniform");expect_identical(st$settings$kernel.admission$kernel,"uniform")
    if(!pooled)expect_identical(st$pilot.bws$ckertype,"gaussian")
    row<-function(e){w<-as.numeric(abs(x$x-e$x)<1.1);if(pooled)w<-w*ifelse(x$group==e$group,.8,.1)
      Q<-outer(x$x-e$x,0:degree,`^`);Z<-Q*sqrt(w);norm<-sqrt(colSums(Z^2));v<-svd(sweep(Z,2L,norm,"/"))
      drop((c(1,rep(0,degree))/norm)%*%v$v%*%(t(v$u)/v$d))*sqrt(w)}
    L<-do.call(rbind,lapply(seq_len(nrow(grid)),function(j)row(grid[j,,drop=FALSE])))
    value<-.np_bias_regression_rows(st,grid,data.frame(row=seq_len(nrow(grid)),target="mean"))
    expect_identical(dim(value$original),dim(L));expect_equal(as.numeric(value$original),as.numeric(L),tolerance=1e-10);expect_equal(value$original_center,drop(L%*%y),tolerance=1e-10)
    expect_equal(value$original_center,fitted(npreg(bw,txdat=x,tydat=y,exdat=grid,se=FALSE)),tolerance=1e-10)
    if(pooled){map<-data.frame(row=seq_len(nrow(grid)),target="c2");got<-.np_bias_regression_rows(st,grid,map)
      expected<-do.call(rbind,lapply(seq_len(nrow(grid)),function(j){e<-grid[j,,drop=FALSE];ref<-e;ref$group<-factor("a",levels=levels(x$group));row(e)-row(ref)}))
      expect_equal(got$original,expected,tolerance=1e-10);expect_true(any(got$se==0))}
    expect_error(.np_bias_regression_rows(st,grid,data.frame(row=seq_len(nrow(grid)),target="d1")),"continuous derivatives remain unqualified")
    expect_error(.np_bias_inference(st,grid,targets="d1",B=19L,law="mammen",calibration="basic-local-envelope",resolution.warn=FALSE),"continuous derivatives remain unqualified")
    expect_identical(bw,held)
  }
})

test_that("public uniform means retain defaults and reject unqualified derivatives", {
  skip_on_cran()
  if(!spawn_mpi_slaves(n=1L))skip("MPI workers unavailable")
  on.exit(close_mpi_slaves(force=TRUE),add=TRUE)
  old<-options(np.messages=FALSE,np.plot.progress=FALSE);on.exit(options(old),add=TRUE)
  for(pooled in c(FALSE,TRUE)) {
    x<-data.frame(x=if(pooled)rep(seq(-.9,.9,length.out=20L),3L)else seq(-.9,.9,length.out=60L))
    if(pooled)x$group<-factor(rep(letters[1:3],each=20L))
    y<-sin(3*x$x)+.15*cos(seq_len(60L))+if(pooled).2*as.integer(x$group)else 0
    bw<-npregbw(xdat=x,ydat=y,bws=c(1.1,if(pooled).2),regtype="lc",ckertype="uniform",bandwidth.compute=FALSE)
    held<-bw;set.seed(112921L);seed<-.Random.seed
    run<-function(control=NULL,derivatives=FALSE){args<-list(x=bw,xdat=x,ydat=y,center="bias-corrected",errors="bootstrap",B=199L,neval=3L,alpha=.1,band="all",gradients=derivatives,output="data");if(!missing(control))args<-c(args,list(boot.control=control));capture_plot_provider(do.call(plot,args))}
    automatic<-run();expect_error(run(NULL),"boot.control must be created by np_boot_control",fixed=TRUE);explicit<-run(np_boot_control(calibration="basic-local-envelope"));expect_identical(automatic$result$panels,explicit$result$panels)
    for(control in list(np_boot_control(),np_boot_control(calibration="default")))expect_identical(run(control)$provider,automatic$provider)
    rad<-run(np_boot_control(wild="rademacher"));tail<-run(np_boot_control(wild="rademacher",calibration="tail-envelope"));expect_identical(rad$result$panels,tail$result$panels)
    expect_identical(rad$result$inference[[1L]]$calibration,"tail-envelope");expect_identical(rad$result$inference[[1L]]$resampling_law,"rademacher")
    expect_error(run(derivatives=TRUE),"continuous derivatives remain unqualified")
    expect_error(run(derivatives=TRUE),"construction='legacy'",fixed=TRUE)
    expect_identical(bw,held);expect_identical(.Random.seed,seed)
  }
})

test_that("fixed rows with saved precision failures preserve live explicit scope", {
  b <- list(type="fixed", regtype="lp", basis="glp", degree=2L, ncon=1L,
            ckertype="epanechnikov", ckerorder=2L, ckerbound="none",
            bernstein.basis=FALSE)
  resolve <- function(b, calibration, gradients=FALSE)
    .np_bias_regression_calibration(b, calibration, "bootstrap", "simultaneous",
                                    "inid", "mammen", gradients)
  for (g in c(FALSE, TRUE)) {
    expect_error(resolve(b, "default", g), "explicit basic-local-envelope")
    expect_identical(resolve(b, "basic-local-envelope", g), "basic-local-envelope")
  }
  b$bernstein.basis <- TRUE
  expect_identical(resolve(b, "default"), "basic-local-envelope")
  b$ckertype <- "beta"
  b$ckerbound <- "fixed"
  expect_identical(resolve(b, "default"), "basic-local-envelope")
  expect_error(resolve(b, "absolute", TRUE), "explicit basic-local-envelope")
  expect_identical(resolve(b, "basic-local-envelope", TRUE), "basic-local-envelope")
  b$degree <- 0L
  expect_identical(resolve(b, "default", TRUE), "basic-local-envelope")
})
