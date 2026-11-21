test_that("partial regression gradients retain undefined panels for every requested band", {
  skip_on_cran()
  opt <- options(np.messages=FALSE,np.plot.progress=FALSE,np.tree=FALSE)
  on.exit(options(opt),add=TRUE)
  local.fit <- function(expr) {
    expr <- substitute(expr)
    if ("npRmpi" %in% loadedNamespaces()) {
      local <- get(".npRmpi_with_local_regression",asNamespace("npRmpi"))
      return(eval(substitute(LOCAL(EXPR),list(LOCAL=local,EXPR=expr)),parent.frame()))
    }
    eval(expr,parent.frame())
  }
  x <- data.frame(x=seq(-.97,.96,length.out=96L),x2=.9*sin(seq_len(96L)*sqrt(2)))
  y <- sin(3*x$x)+.1*cos(seq_len(96L))
  set.seed(1127701);seed <- .Random.seed
  for(type in c("fixed","generalized_nn","adaptive_nn"))for(degree in list(c(0L,1L),c(1L,0L))) {
    bw <- local.fit(npregbw(xdat=x,ydat=y,bws=if(type=="fixed")c(.5,.5)else c(45,45),
      regtype="lp",degree=degree,basis="glp",bernstein.basis=TRUE,bwtype=type,ckertype="gaussian",bandwidth.compute=FALSE))
    fit <- local.fit(suppressWarnings(npreg(bw,txdat=x,tydat=y,gradients=TRUE)));held <- fit
    for(neval in c(1L,3L)) {
      z <- suppressWarnings(plot(fit,errors="bootstrap",center="bias-corrected",bootstrap="wild",
        B=63L,random.seed=1127702L,gradients=TRUE,neval=neval,alpha=.2,band="all",output="data"))
      expect_s3_class(z,"np_plot_data")
      missing <- z$panels[[which(degree==0L)]];available <- z$panels[[which(degree==1L)]]
      expect_equal(nrow(missing$values),neval)
      expect_true(all(missing$values$status=="undefined"))
      expect_true(all(is.na(missing$values$estimate)))
      expect_true(all(is.na(missing$values$bias_corrected)))
      expect_equal(nrow(missing$intervals),3L*neval)
      expect_setequal(missing$intervals$band,c("pointwise","simultaneous","bonferroni"))
      expect_true(all(missing$intervals$status=="undefined"))
      expect_true(all(is.na(missing$intervals$lower)))
      expect_true(all(is.na(missing$intervals$upper)))
      expect_true(all(available$values$status=="available"))
      expect_true(all(is.finite(available$intervals$lower)))
      expect_true(all(is.finite(available$intervals$upper)))
    }
    expect_identical(fit,held);expect_identical(.Random.seed,seed)
  }
})
