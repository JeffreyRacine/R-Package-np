a2_capture <- function(expr) {
  warning.state <- new.env(hash = FALSE, parent = emptyenv())
  warning.state$messages <- character()
  value <- withCallingHandlers(tryCatch(expr, error = identity), warning = function(w) {
    warning.state$messages <- c(warning.state$messages, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  list(value = value, warnings = warning.state$messages[grepl("all computed kernel weights", warning.state$messages)])
}

a2_fixture <- function() {
  z <- data.frame(z = c(seq(-1, -.5, length.out = 32), seq(.5, 1, length.out = 32)))
  x <- data.frame(x = sin(seq_len(64)))
  list(z=z, x=x, y=sin(z$z)+.3*x$x+.1*cos(seq_len(64)),
       ez=data.frame(z=c(-.75, 0, .75)), ex=data.frame(x=c(.1, .2, .3)))
}

a2_expect_plot <- function(expr) {
  got <- a2_capture(expr)
  expect_false(inherits(got$value, "error"))
  expect_length(got$warnings, 1L)
  expect_match(got$warnings, "grid/row")
  expect_null(attr(got$value, ".np.empty.rows", exact=TRUE))
  invisible(got$value)
}

test_that("regression plots collect undefined rows over complete public grids", {
  old <- options(np.messages=FALSE, np.tree=FALSE); on.exit(options(old), add=TRUE)
  d <- a2_fixture()
  for (spec in list(list(regtype="ll"), list(regtype="lp", degree=3L, bernstein.basis=TRUE))) {
    b <- do.call(npregbw, c(list(xdat=d$z, ydat=d$y, bws=.12,
      bandwidth.compute=FALSE, ckertype="epanechnikov"), spec))
    a2_expect_plot(plot(b, xdat=d$z, ydat=d$y, neval=9,
      plot.behavior="data", plot.errors.method="none"))
    for (method in c("none", "asymptotic")) {
      a2_expect_plot(plot(b, xdat=d$z, ydat=d$y, neval=9,
        plot.behavior="data", plot.errors.method=method, gradients=TRUE))
    }
    helper <- a2_capture(.np_plot_regression_eval(b, d$z, d$y, d$ez, gradients=TRUE))
    expect_length(helper$warnings, 1L)
    expect_true(is.na(helper$value$mean[2]))
    expect_true(is.na(helper$value$grad[2,1]))
    ref <- .np_plot_regression_eval(b, d$z, d$y, d$ez[c(1,3),], gradients=TRUE)
    expect_equal(helper$value$mean[c(1,3)], ref$mean, tolerance=1e-12)
    expect_equal(helper$value$grad[c(1,3),,drop=FALSE], ref$grad, tolerance=1e-12)
    expect_error(.np_regression_direct(b, d$z, d$y, d$ez), "LP solve")
    expect_error(plot(b, .np.empty.report=function(...)NULL), "np.empty.report|unknown|unsupported|unrecognized")
  }
  z <- expand.grid(z1=d$z$z[c(1:8,33:40)], z2=d$z$z[c(1:8,33:40)])
  y <- sin(z$z1)+z$z2^2
  b <- npregbw(xdat=z, ydat=y, bws=c(.12,.12), bandwidth.compute=FALSE,
    regtype="lp", degree=c(1L,3L), bernstein.basis=TRUE, ckertype="epanechnikov")
  got <- a2_capture(plot(b, xdat=z, ydat=y, neval=9, perspective=FALSE,
    plot.behavior="data", plot.errors.method="none"))
  expect_false(inherits(got$value,"error")); expect_length(got$warnings,1L)
  expect_match(got$warnings,"18 external evaluation row")
})

test_that("single-index and partially linear plot owners keep strict internals", {
  old <- options(np.messages=FALSE, np.tree=FALSE); on.exit(options(old), add=TRUE)
  d <- a2_fixture()
  b <- suppressWarnings(npindexbw(xdat=d$z, ydat=d$y, bws=c(1,.12),
    bandwidth.compute=FALSE, regtype="lp", degree=3L, ckertype="epanechnikov"))
  for (method in c("none","asymptotic"))
    a2_expect_plot(plot(b, xdat=d$z, ydat=d$y, neval=9,
      plot.behavior="data", plot.errors.method=method, gradients=TRUE))
  for (type in c("generalized_nn","adaptive_nn")) {
    nn <- suppressWarnings(npindexbw(xdat=d$z, ydat=d$y, bws=c(1,16),
      bwtype=type, bandwidth.compute=FALSE, regtype="lp", degree=3L))
    idx <- data.frame(index=d$z$z)
    result <- .np_plot_singleindex_local_eval(nn, idx, idx[c(10,50),,drop=FALSE], d$y)
    ref <- .np_indexhat_exact(nn, idx, idx[c(10,50),,drop=FALSE], d$y, output="apply")
    expect_equal(result$mean, as.vector(ref), tolerance=1e-12)
    expect_null(attr(result,".np.empty.rows",exact=TRUE))
  }
  bp <- npplregbw(xdat=d$x, zdat=d$z, ydat=d$y, bws=matrix(.12,2,1),
    bandwidth.compute=FALSE, regtype="lp", degree=3L, ckertype="epanechnikov")
  for (method in c("none","asymptotic"))
    a2_expect_plot(plot(bp, xdat=d$x, zdat=d$z, ydat=d$y, neval=9, perspective=FALSE,
      plot.behavior="data", plot.errors.method=method))
  h <- a2_capture(.np_plot_plreg_local_fit(bp,d$x,d$y,d$z,d$ex,d$ez))
  expect_length(h$warnings,1L); expect_true(is.na(h$value$mean[2]))
  ref <- .np_plot_plreg_local_fit(bp,d$x,d$y,d$z,d$ex[c(1,3),,drop=FALSE],d$ez[c(1,3),,drop=FALSE])
  expect_equal(h$value$mean[c(1,3)],ref$mean,tolerance=1e-12)
})

test_that("vector least-squares quantile prediction and plots publish once", {
  old <- options(np.messages=FALSE,np.tree=FALSE); on.exit(options(old),add=TRUE)
  d <- a2_fixture()
  b <- nplsqregbw(xdat=d$z,ydat=d$y,tau=c(.25,.75),scale=rep(.3,64),bws=.12,
    bandwidth.compute=FALSE,regtype="lp",degree=3L,ckertype="epanechnikov")
  f <- nplsqreg(b)
  for (se.fit in c(FALSE,TRUE)) {
    got <- a2_capture(predict(f,exdat=d$ez,se.fit=se.fit))
    expect_false(inherits(got$value,"error")); expect_length(got$warnings,1L)
    val <- if(se.fit)got$value$fit else got$value
    expect_true(all(is.na(val[2,]))); expect_true(all(is.finite(val[c(1,3),])))
    expect_null(attr(got$value,".np.empty.rows",exact=TRUE))
    ref <- predict(f,exdat=d$ez[c(1,3),,drop=FALSE],se.fit=se.fit)
    expect_equal(val[c(1,3),,drop=FALSE],if(se.fit)ref$fit else ref,tolerance=1e-12)
  }
  for (method in c("none","asymptotic"))
    a2_expect_plot(plot(f,neval=9,plot.behavior="data",plot.errors.method=method))
  options(warn=2)
  expect_error(predict(f,exdat=d$ez),"all computed kernel weights")
  expect_error(plot(f,neval=9,plot.behavior="data",plot.errors.method="none"),"all computed kernel weights")
  options(warn=0)
  expect_true(all(is.finite(predict(f,exdat=d$ez[c(1,3),,drop=FALSE]))))
})


test_that("a required bootstrap failure publishes no partial plot warning", {
  # Inject at the required public bootstrap stage after the point grid has run.
  # This tests publication ordering without changing a numerical failure policy.
  local_mocked_bindings(
    compute.bootstrap.errors.rbandwidth = function(...) stop("required bootstrap witness"),
    .package = "npRmpi"
  )
  old <- options(np.messages=FALSE,np.tree=FALSE); on.exit(options(old),add=TRUE)
  d <- a2_fixture()
  b <- npregbw(xdat=d$z,ydat=d$y,bws=.12,bandwidth.compute=FALSE,
    regtype="lp",degree=3L,ckertype="epanechnikov")
  got <- a2_capture(plot(b,xdat=d$z,ydat=d$y,neval=9,plot.behavior="data",
    plot.errors.method="bootstrap",plot.errors.boot.num=3))
  expect_s3_class(got$value,"error")
  expect_match(conditionMessage(got$value),"required bootstrap witness")
  expect_length(got$warnings,0L)
  a2_expect_plot(plot(b,xdat=d$z,ydat=d$y,neval=9,plot.behavior="data",
    plot.errors.method="none"))
})
