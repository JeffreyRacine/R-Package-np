
test_that("fixed uniform degree-zero main rows preserve physical support ties", {
  skip_on_cran()
  opt <- options(np.messages=FALSE,np.plot.progress=FALSE,np.tree=FALSE)
  on.exit(options(opt),add=TRUE)
  set.seed(4202)
  x <- data.frame(x=c(18,67,sample(18:67,148L,replace=TRUE)))
  y <- log(x$x)+rnorm(150L,sd=.1)
  seed <- .Random.seed
  for(bound in c("none","fixed")) for(reg in c("lc","lp")) {
    bw <- suppressWarnings(npregbw(xdat=x,ydat=y,bws=3,bwtype="fixed",
      ckertype="uniform",ckerbound=bound,ckerlb=if(bound=="fixed")18 else NULL,
      ckerub=if(bound=="fixed")67 else NULL,regtype=reg,degree=0L,bandwidth.compute=FALSE))
    fit <- npreg(bw,txdat=x,tydat=y); held <- fit
    value <- suppressWarnings(capture_plot_provider(plot(fit,center="bias-corrected",
      errors="bootstrap",bootstrap="inid",neval=50L,B=31L,
      random.seed=1120501L,output="data")))
    inf <- value$provider[[1L]]$inference
    grid <- inf$grid
    # Independent physical-coordinate LC weights and query-common bounded
    # mass. The strict open support convention is unchanged.
    cdf <- function(z).5+.5*pmax(-1,pmin(1,z))
    expected <- vapply(grid$x,function(g) {
      h <- 3
      at <- g
      mass <- if(bound=="fixed")cdf((67-at)/h)-cdf((18-at)/h) else 1
      w <- .5*(abs((x$x-g)/h)<1)/h/mass
      sum(w*y)/sum(w)
    },0)
    expect_equal(inf$estimate,expected,tolerance=1e-12)
    expect_equal(inf$estimate,npreg(bw,txdat=x,tydat=y,exdat=grid)$mean,tolerance=1e-12)
    expect_identical(fit,held);expect_identical(.Random.seed,seed)
  }
})
