test_that("physical bandwidth extraction preserves scaled profile models", {
  old <- options(np.messages=FALSE,np.categorical.compress=TRUE)
  on.exit(options(old),add=TRUE)
  i <- seq_len(30L)
  x <- data.frame(x=sin(i*.73))
  z <- data.frame(u=factor(rep(letters[1:3],10)),o=ordered(rep(1:2,15)))
  y <- cos(i*.27)+x$x*as.integer(z$u)
  for (family in c("npreg","npscoef","npplreg")) {
    a <- if(family=="npreg") list(xdat=z,ydat=y) else list(xdat=x,zdat=z,ydat=y)
    raw <- if(family=="npplreg") matrix(c(.3,.4),2,2,byrow=TRUE) else c(.3,.4)
    make <- get(paste0(family,"bw"))
    seed <- do.call(make,c(a,list(bws=raw,bwscaling=TRUE,bandwidth.compute=FALSE)))
    factor <- if(family=="npplreg") seed$bw$yzbw$ncatfac else seed$ncatfac
    scaled <- do.call(make,c(a,list(bws=raw/factor,bwscaling=TRUE,bandwidth.compute=FALSE)))
    physical <- do.call(make,c(a,list(bws=raw,bwscaling=FALSE,bandwidth.compute=FALSE)))
    fit <- function(b) if(family=="npscoef") npscoef(b,iterate=FALSE) else get(family)(b)
    fs <- fit(scaled); fp <- fit(physical)
    expect_equal(fitted(fs),fitted(fp),tolerance=2e-11)
    if(family!="npreg") expect_equal(coef(fs),coef(fp),tolerance=2e-11)
    if(family=="npreg") {
      codes <- .np_cat_profile_code_matrix(z)
      W <- .np_regression_cat_profile_kernel_matrix(codes,codes,z,scaled)
      # Independent finite-support weights: unordered Aitchison-Aitken and
      # ordered LR, the regression-family defaults.
      oracle <- ifelse(outer(z$u,z$u,"=="),.7,.15) *
        .4^abs(outer(as.integer(z$o),as.integer(z$o),"-"))
      expect_equal(W,oracle,tolerance=1e-14)
      expect_equal(as.double(fitted(fs)),as.double(W%*%y/rowSums(W)),tolerance=2e-12)
    }
  }
})

test_that("scaled hats and bootstrap reuse equal their physical estimator", {
  old <- options(np.messages=FALSE)
  on.exit(options(old),add=TRUE)
  i <- seq_len(30L)
  x <- data.frame(x=sin(i*.73)+1.2,g=factor(rep(letters[1:3],10)))
  y <- cos(i*.27)+.2*x$x
  ex <- x[c(3,7,11,18),]
  for(reg in c("lc","ll","lp")) {
    a <- list(xdat=x,ydat=y,bws=c(1.3,.4),bwscaling=TRUE,
              bandwidth.compute=FALSE,regtype=reg)
    if(reg=="lp") a$degree <- 2L
    sc <- do.call(npregbw,a)
    a$bws <- unlist(sc$bandwidth); a$bwscaling <- FALSE
    ph <- do.call(npregbw,a)
    for(s in 0:1) {
      hs <- npreghat(sc,txdat=x,exdat=ex,s=s,output="matrix")
      hp <- npreghat(ph,txdat=x,exdat=ex,s=s,output="matrix")
      expect_equal(as.double(hs),as.double(hp),tolerance=2e-11)
      ap <- npreghat(sc,txdat=x,exdat=ex,s=s,y=y,output="apply")
      expect_equal(as.double(ap),as.double(hp%*%y),tolerance=2e-11)
      f <- npreg(ph,txdat=x,tydat=y,exdat=ex,gradients=TRUE)
      ref <- if(s==0L) fitted(f) else gradients(f)[,1L]
      expect_equal(as.double(hs%*%y),as.double(ref),tolerance=2e-10)
    }
    run <- function(b) {
      set.seed(819)
      p <- plot(b,xdat=x,ydat=y,plot.behavior="data",
                plot.errors.method="bootstrap",plot.errors.boot.num=9L,
                plot.errors.boot.method="wild",neval=7L)
      list(values=lapply(p,function(z)c(z$mean,z$merr)),rng=.Random.seed)
    }
    ss <- run(sc); pp <- run(ph)
    expect_equal(ss$values,pp$values,tolerance=2e-10)
    expect_identical(ss$rng,pp$rng)
  }
})
