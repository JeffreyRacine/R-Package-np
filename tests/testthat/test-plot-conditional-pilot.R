cf75_get <- function(name) getFromNamespace(name, "npRmpi")
cf75_local <- function(expr) {
  if ("npRmpi" == "npRmpi") cf75_get(".npRmpi_with_local_regression")(force(expr))
  else force(expr)
}
cf75_fixture <- function(cdf=FALSE, regtype="ll", kernel="gaussian", ordered=FALSE) {
  i <- seq_len(36L)
  x <- data.frame(f=factor(rep(c("a","b","c"),12), ordered=ordered),
                  x=sin(i*.73))
  y <- data.frame(y=cos(i*.27)+.2*x$x)
  constructor <- cf75_get(if(cdf) "npcdistbw" else "npcdensbw")
  args <- list(xdat=x, ydat=y, bws=c(.6,.2,.7),
    cxkertype=kernel, cykertype=kernel, regtype=regtype,
    bandwidth.compute=FALSE)
  if(regtype=="lp") args$degree <- 2L
  bw <- cf75_local(do.call(constructor,args))
  list(x=x,y=y,bw=bw,ex=x[c(1,2,3,7,8,9),],
       ey=data.frame(y=rep(.1,6L)))
}

test_that("smooth perturbations use estimation units without new RNG consumption", {
  random <- cf75_get(".np_plot_kernel_random")
  old.epan <- cf75_get(".np_plot_repan2")
  set.seed(75); old <- old.epan(200L); seed <- .Random.seed
  set.seed(75); new <- random(200L,"epanechnikov",2L)
  expect_identical(new,sqrt(5)*old)
  expect_identical(.Random.seed,seed)
  for (kernel in c("gaussian","uniform")) {
    set.seed(75)
    expected <- if(kernel=="gaussian") rnorm(200L) else runif(200L,-1,1)
    seed <- .Random.seed
    set.seed(75)
    expect_identical(random(200L,kernel,2L),expected)
    expect_identical(.Random.seed,seed)
  }
  expect_error(random(2L,"gaussian",4L),"signed")
  expect_error(random(2L,"epanechnikov",4L),"signed")
  expect_length(random(0L,"epanechnikov",2L),0L)
})

test_that("the conditional pilot has its own coherent truth and paired contrast", {
  prepare <- cf75_get(".np_plot_conditional_pilot_prepare")
  reference <- cf75_get(".np_plot_conditional_pilot_reference")
  oversmooth <- cf75_get(".np_plot_oversmooth_conditional_bws")
  for (cdf in c(FALSE,TRUE)) {
    f <- cf75_fixture(cdf)
    f$ex$x <- seq(-.4,.4,length.out=6L)
    g <- oversmooth(f$bw,cdf)
    pilot <- prepare(f$x,f$y,f$bw,cdf)
    oracle <- function(ex) vapply(seq_len(nrow(ex)), function(j) {
      w <- ifelse(f$x$f==ex$f[j],1-g$bandwidth$x[1],g$bandwidth$x[1]/2) *
        dnorm((ex$x[j]-f$x$x)/g$bandwidth$x[2])
      v <- if(cdf) pnorm((f$ey$y[j]-f$y$y)/g$bandwidth$y[1]) else
        dnorm((f$ey$y[j]-f$y$y)/g$bandwidth$y[1])/g$bandwidth$y[1]
      sum(w*v)/sum(w)
    }, numeric(1))
    expect_equal(reference(pilot,f$ex,f$ey),oracle(f$ex),tolerance=1e-13)
    lower <- f$ex; lower$f <- factor("a",levels=levels(f$x$f))
    expect_equal(reference(pilot,f$ex,f$ey,1L),oracle(f$ex)-oracle(lower),
                 tolerance=1e-13)
  }
})

test_that("categorical pilot normalization preserves ordered support and limits", {
  side <- cf75_get(".np_plot_pilot_prepare_side")
  raw <- cf75_get(".np_plot_pilot_cat_raw")
  draw <- cf75_get(".np_plot_pilot_draw_side")
  x <- data.frame(o=ordered(rep(c("1","3","8"),4L),levels=c("1","3","8")))
  for(kernel in c("wangvanryzin","liracine","nliracine","racineliyan"))
    for(lambda in c(0,.3,1)) {
      p <- side(x,lambda,FALSE,"gaussian",2L,NULL,NULL,"aitchisonaitken",kernel,TRUE)
      s <- p$specs[[1L]]
      Q <- raw(s,s$support,s$support)
      Q <- sweep(Q,2L,s$mass,"/")
      expect_equal(colSums(Q),rep(1,3),tolerance=1e-14)
      expect_true(all(Q>=0))
      if(lambda==0) expect_equal(Q,diag(3))
      set.seed(75); result <- draw(p,rep(1:3,each=40L))$o
      expect_identical(levels(result),levels(x$o))
      expect_true(is.ordered(result))
      if(lambda==0) expect_identical(as.character(result),rep(c("1","3","8"),each=40L))
    }
})

test_that("conditional smooth refits honor LC LL LP and physical gradient targets", {
  boot <- cf75_get(".np_plot_conditional_pilot_boot")
  eval <- cf75_get(".np_plot_conditional_eval")
  for(cdf in c(FALSE,TRUE)) for(regtype in c("lc","ll","lp"))
    for(ordered in c(FALSE,TRUE)) {
      f <- cf75_fixture(cdf,regtype,ordered=ordered)
      for(target in list(NULL,1L,2L)) {
        args <- list(xdat=f$x,ydat=f$y,exdat=f$ex,eydat=f$ey,bws=f$bw,cdf=cdf,
          plot.errors.boot.method="inid",plot.errors.boot.blocklen=NULL,
          plot.errors.boot.num=2L,progress.label=NULL,gradient.index=target)
        set.seed(751L)
        result <- cf75_local(do.call(boot,args))
        expected <- cf75_local(eval(f$bw,f$x,f$y,f$ex,f$ey,cdf=cdf,
          gradients=!is.null(target),se=FALSE,gradient.target=target,
          lp.first.se.demand=FALSE,cat.se.demand=FALSE))
        value <- if(is.null(target)) expected[[if(cdf) "condist" else "condens"]] else
          expected$congrad[,target]
        expect_identical(as.vector(result$t0),as.vector(value))
        expect_identical(dim(result$t),c(2L,6L))
        expect_true(all(is.finite(result$t)))
        expect_true(all(is.finite(result$center)))
      }
    }
})

test_that("finite-bound continuous pilot draws and reference share one law", {
  side <- cf75_get(".np_plot_pilot_prepare_side")
  draw <- cf75_get(".np_plot_pilot_draw_side")
  pc <- cf75_get(".np_plot_pilot_continuous")
  pq <- cf75_get(".np_plot_pilot_quantile")
  for(kernel in c("gaussian","uniform","epanechnikov")) {
    p <- seq(.001,.999,length.out=21L)
    expect_equal(pc(pq(p,kernel),kernel,cdf=TRUE),p,tolerance=1e-13)
    s <- side(data.frame(x=c(.01,.5,.99)),.7,TRUE,kernel,2L,0,1,NULL,NULL,FALSE)
    set.seed(752L)
    z <- draw(s,rep(1:3,each=40L))$x
    expect_true(all(z>=0 & z<=1))
    for(donor in 1:3) {
      radius <- if(kernel=="gaussian") Inf else if(kernel=="uniform") .7 else .7*sqrt(5)
      mass <- integrate(function(v) pc((v-s$codes[donor,1])/.7,kernel)/
        (.7*s$specs[[1L]]$mass[donor]),max(0,s$codes[donor,1]-radius),
        min(1,s$codes[donor,1]+radius),rel.tol=1e-12)$value
      expect_equal(mass,1,tolerance=1e-11)
    }
  }
})
test_that("higher pilot derivatives use physical coordinates and selected order", {
  f <- cf75_fixture(FALSE,"lp")
  f$x$z <- cos(seq_len(nrow(f$x))*.17)
  f$ex <- f$x[c(1,2,3,7,8,9),]
  f$ex$x <- seq(-.4,.4,length.out=6L)
  f$bw <- cf75_local(cf75_get("npcdensbw")(xdat=f$x,ydat=f$y,bws=c(.6,.2,.7,.8),
    regtype="lp",degree=c(2L,2L),bandwidth.compute=FALSE))
  g <- cf75_get(".np_plot_oversmooth_conditional_bws")(f$bw,FALSE)
  pilot <- cf75_get(".np_plot_conditional_pilot_prepare")(f$x,f$y,f$bw,FALSE)
  expected <- vapply(seq_len(nrow(f$ex)),function(i) {
    z <- (f$ex$x[i]-f$x$x)/g$bandwidth$x[2]
    a <- ifelse(f$x$f==f$ex$f[i],1-g$bandwidth$x[1],g$bandwidth$x[1]/2) *
      dnorm(z)*dnorm((f$ex$z[i]-f$x$z)/g$bandwidth$x[3])
    a1 <- -z*a/g$bandwidth$x[2]
    a2 <- (z*z-1)*a/g$bandwidth$x[2]^2
    b <- dnorm((f$ey$y[i]-f$y$y)/g$bandwidth$y[1])/g$bandwidth$y[1]
    D <- sum(a); D1 <- sum(a1); D2 <- sum(a2)
    N <- sum(a*b); N1 <- sum(a1*b); N2 <- sum(a2*b)
    N2/D-2*D1*N1/D^2-D2*N/D^2+2*D1^2*N/D^3
  },numeric(1))
  reference <- cf75_get(".np_plot_conditional_pilot_reference")
  expect_equal(reference(pilot,f$ex,f$ey,2L,2L),expected,tolerance=1e-12)
  result <- cf75_local(cf75_get(".np_plot_conditional_pilot_boot")(
    f$x,f$y,f$ex,f$ey,f$bw,FALSE,"inid",NULL,2L,NULL,
    gradient.index=2L,gradient.order=c(2L,1L)))
  expect_equal(result$center,expected,tolerance=1e-12)
  expect_true(all(is.finite(result$t)))
})

test_that("mixed response pilot atoms and scalar kernels define one PDF CDF law", {
  for(cdf in c(FALSE,TRUE)) {
    f <- cf75_fixture(cdf,"lc")
    f$y$v <- factor(rep(c("left","right"),18L))
    f$ey <- f$y[c(1,2,3,7,8,9),]
    f$ey$y <- .1
    f$bw <- cf75_local(cf75_get(if(cdf) "npcdistbw" else "npcdensbw")(
      xdat=f$x,ydat=f$y,bws=c(.6,.15,.2,.7),regtype="lc",
      cykertype="epanechnikov",bandwidth.compute=FALSE))
    pilot <- cf75_get(".np_plot_conditional_pilot_prepare")(f$x,f$y,f$bw,cdf)
    g <- cf75_get(".np_plot_oversmooth_conditional_bws")(f$bw,cdf)
    expected <- vapply(seq_len(nrow(f$ex)),function(i) {
      a <- ifelse(f$x$f==f$ex$f[i],1-g$bandwidth$x[1],g$bandwidth$x[1]/2) *
        dnorm((f$ex$x[i]-f$x$x)/g$bandwidth$x[2])
      z <- (f$ey$y[i]-f$y$y)/g$bandwidth$y[1]/sqrt(5)
      u <- pmax(-1,pmin(1,z))
      b <- if(cdf) .5+.75*u-.25*u^3 else
        ifelse(abs(z)<1,3*(1-z^2)/(4*sqrt(5)*g$bandwidth$y[1]),0)
      atom <- if(cdf) as.integer(f$y$v)<=as.integer(f$ey$v[i]) else f$y$v==f$ey$v[i]
      sum(a*b*atom)/sum(a)
    },numeric(1))
    expect_equal(cf75_get(".np_plot_conditional_pilot_reference")(pilot,f$ex,f$ey),
                 expected,tolerance=1e-12)
  }
})
test_that("fixed and geometric pilots retain paired targets and compact kernel refits", {
  for(cdf in c(FALSE,TRUE)) for(kernel in c("epanechnikov","uniform"))
    for(method in c("fixed","geom")) {
      known.warnings <- 0L
      out <- withCallingHandlers({
        f <- cf75_fixture(cdf,"ll",kernel,ordered=TRUE)
        set.seed(753L)
        cf75_local(cf75_get(".np_plot_conditional_pilot_boot")(
          f$x,f$y,f$ex,f$ey,f$bw,cdf,method,3L,3L,NULL,gradient.index=1L))
      }, warning=function(w) {
        if(identical(conditionMessage(w),
          "[npRmpi] ignoring kernel order specified with uniform kernel type")) {
          known.warnings <<- known.warnings+1L
          invokeRestart("muffleWarning")
        }
      })
      expect_identical(known.warnings>0L,kernel=="uniform")
      expect_identical(dim(out$t),c(3L,6L))
      expect_true(all(is.finite(out$t)))
      expect_true(all(is.finite(out$center)))
    }
})

test_that("conditional pilot handles categorical-only X and bounded selected fits", {
  for(cdf in c(FALSE,TRUE)) {
    f <- cf75_fixture(cdf,"lc")
    f$x <- f$x["f"]; f$ex <- f$ex["f"]
    constructor <- cf75_get(if(cdf) "npcdistbw" else "npcdensbw")
    f$bw <- cf75_local(constructor(xdat=f$x,ydat=f$y,bws=c(.6,.2),
      bandwidth.compute=FALSE))
    out <- cf75_local(cf75_get(".np_plot_conditional_pilot_boot")(
      f$x,f$y,f$ex,f$ey,f$bw,cdf,"inid",NULL,2L,NULL,gradient.index=1L))
    expect_true(all(is.finite(out$t)))
    expect_true(all(is.finite(out$center)))
  }
  f <- cf75_fixture(FALSE,"ll")
  f$bw <- cf75_local(cf75_get("npcdensbw")(xdat=f$x,ydat=f$y,
    bws=c(.6,.2,.7),regtype="ll",cxkerbound="fixed",cxkerlb=c(NA,-1.1),
    cxkerub=c(NA,1.1),cykerbound="fixed",cykerlb=-1.5,cykerub=1.5,
    bandwidth.compute=FALSE))
  out <- cf75_local(cf75_get(".np_plot_conditional_pilot_boot")(
    f$x,f$y,f$ex,f$ey,f$bw,FALSE,"inid",NULL,3L,NULL,gradient.index=1L))
  expect_true(all(is.finite(out$t)))
  expect_true(all(is.finite(out$center)))
})
