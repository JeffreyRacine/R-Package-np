test_that("quantile plot replay retains controls and explicit precedence", {
  pkg <- "np"
  old <- options(np.messages = FALSE); on.exit(options(old), add = TRUE)
  d <- data.frame(x = seq(-1, 1, length.out = 15L))
  d$y <- d$x + .2*cos(seq_len(nrow(d)))
  b <- npcdistbw(y ~ x, data = d, bws = c(.4, .6), bandwidth.compute = FALSE)
  ctl <- list(tol = 1e-9, small = 1e-10, itmax = 80L)
  for (tau in list(.4, c(.3, .7))) for (grad in c(FALSE, TRUE)) {
    fit <- do.call(npqreg, c(list(bws = b, tau = tau), ctl))
    out <- plot(fit, output = "data", errors = "none", perspective = FALSE,
                gradients = grad, neval = 3L)
    expect_length(out, 1L)
    panel <- out[[1L]]
    ref <- do.call(npqreg, c(list(bws = b, exdat = panel$xeval,
                                  tau = tau, gradients = grad), ctl))
    expect_identical(panel$fit.controls, ctl)
    expect_equal(as.vector(panel$quantile), as.vector(ref$quantile), tolerance = 1e-12)
    if (grad) expect_equal(as.vector(panel$quantgrad), as.vector(ref$quantgrad),
                          tolerance = 1e-12)
    over <- plot(fit, output = "data", errors = "none", perspective = FALSE,
                 gradients = grad, neval = 3L, tol = .01)[[1L]]
    expect_identical(over$fit.controls, modifyList(ctl, list(tol = .01)))
  }
  expect_error(plot(fit, output = "data", tol = NULL), "'tol'")
  expect_error(plot(fit, output = "data", small = NULL), "'small'")
  expect_error(plot(fit, output = "data", itmax = 1L), "failed to converge")
  legacy <- fit; legacy$fit.controls <- NULL
  a <- plot(legacy, output = "data", neval = 3L)[[1L]]
  expect_identical(a$fit.controls,
    list(tol = 1.490116e-04, small = 1.490116e-05, itmax = 10000L))
  expect_error(plot(b, output = "data", quantreg = FALSE, tol = .01),
               "only to quantile plots")
})

test_that("quantile resamples use controls for levels and categorical gradients", {
  pkg <- "np"
  old <- options(np.messages = FALSE); on.exit(options(old), add = TRUE)
  ns <- asNamespace(pkg)
  local <- if (pkg == "npRmpi") get(".npRmpi_with_local_bootstrap", ns) else
    function(expr) force(expr)
  x <- data.frame(x = seq(-1, 1, length.out = 16L),
                  u = factor(rep(c("a", "b"), 8L)))
  y <- x$x + .2*as.integer(x$u) + .15*sin(seq_len(nrow(x)))
  b <- npcdistbw(xdat = x, ydat = y, bws = c(.5, .7, .2),
                 bandwidth.compute = FALSE)
  ex <- x[c(4L, 9L), , drop = FALSE]
  ctl <- list(tol = 1e-9, small = 1e-10, itmax = 80L)
  counts <- cbind(rep(1L,16L), c(2L,0L,rep(1L,14L)),
                   c(rep(1L,14L),0L,2L))
  ev <- get(".np_plot_quantile_eval", ns)
  for (tau in list(.4, c(.3,.7))) for (component in 0:2) {
    helper <- get(paste0(if (pkg=="npRmpi") ".npRmpi_inid_boot_from_quantile_" else
                        ".np_inid_boot_from_quantile_",
                        if (component==0L) "level" else "gradient",
                        if (pkg=="np") "_local" else ""), ns)
    a <- c(list(xdat=x, ydat=y, exdat=ex, bws=b, B=3L, tau=tau, counts=counts), ctl)
    if (component>0L) a$gradient.index <- component
    out <- do.call(helper,a)
    oracle <- function(idx) {
      z <- local(do.call(ev,c(list(bws=b,txdat=x[idx,,drop=FALSE],tydat=y[idx],
              exdat=ex,tau=tau,gradients=component>0L,need.errors=FALSE,
              lp.first.se.demand=FALSE,cat.se.demand=FALSE),ctl)))
      if(component==0L) as.vector(z$quantile) else
        if(length(dim(z$quantgrad))==3L) as.vector(z$quantgrad[,component,]) else
          as.vector(z$quantgrad[,component])
    }
    expect_identical(out$t0,oracle(seq_len(nrow(x))))
    ref <- t(vapply(seq_len(ncol(counts)),function(j)
      oracle(rep.int(seq_len(nrow(x)),counts[,j])),numeric(length(out$t0))))
    expect_identical(out$t,ref)
  }
})

test_that("quantile public bootstrap carries controls across supported block paths", {
  pkg <- "np"
  old <- options(np.messages = FALSE); on.exit(options(old), add = TRUE)
  d <- data.frame(x=seq(-1,1,length.out=16L))
  d$y <- d$x + .2*cos(seq_len(nrow(d)))
  b <- npcdistbw(y~x,data=d,bws=c(.4,.6),bandwidth.compute=FALSE)
  ctl <- list(tol=1e-9,small=1e-10,itmax=80L)
  fit <- do.call(npqreg,c(list(bws=b,tau=c(.3,.7)),ctl))
  for (method in c("inid","fixed","geom")) {
    set.seed(815)
    out <- plot(fit,output="data",errors="bootstrap",bootstrap=method,
                B=3L,plot.errors.boot.blocklen=2L,neval=3L,
                center="estimate",perspective=FALSE)
    expect_length(out,1L)
    expect_identical(out[[1L]]$fit.controls,ctl)
    expect_true(all(is.finite(out[[1L]]$quantile)))
    expect_true(all(is.finite(out[[1L]]$quanterr)))
    expect_error(plot(fit,output="data",errors="bootstrap",bootstrap=method,
      B=3L,neval=3L,center="bias-corrected"), 'not supported')
  }
})
