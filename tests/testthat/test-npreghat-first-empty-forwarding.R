p07_capture <- function(expr) {
  observed <- new.env(parent=emptyenv())
  observed$warnings <- character()
  value <- withCallingHandlers(expr, warning=function(w) {
    observed$warnings <- c(observed$warnings,conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  list(value=value,warnings=observed$warnings)
}

test_that("first LP hats share the external empty-row policy across geometry", {
  old <- options(np.messages=FALSE,np.tree=FALSE)
  on.exit(options(old),add=TRUE)
  x <- data.frame(x=seq(.1,.9,length.out=41),
                  f=factor(rep(c("a","b"),length.out=41),levels=c("a","b","c")))
  ex <- data.frame(x=c(.25,.25,.65),f=factor(c("a","c","b"),levels=levels(x$f)))
  for(kernel in c("epanechnikov","gaussian"))
    for(type in c("fixed","generalized_nn","adaptive_nn"))
      for(degree in c(1L,2L)) {
        y <- if(degree==1L) 2+3*x$x else x$x^2
        bw <- npregbw(xdat=x,ydat=y,bws=c(if(type=="fixed") .18 else 29,0),
          regtype="lp",degree=degree,ckertype=kernel,bwtype=type,
          ckerbound="fixed",ckerlb=0,ckerub=1,bandwidth.compute=FALSE)
        h <- p07_capture(npreghat(bw,txdat=x,exdat=ex,s=1))
        H <- h$value
        expect_length(h$warnings,1L)
        expect_match(h$warnings,"all computed kernel weights are zero")
        expect_true(all(is.na(H[2,])))
        expect_true(all(is.finite(H[c(1,3),])))
        expect_true(is.na(attr(H,"ridge.used")[2]))
        expect_null(attr(H,".np.empty.rows",exact=TRUE))
        truth <- if(degree==1L) c(3,3) else 2*ex$x[c(1,3)]
        expect_equal(drop(H[c(1,3),,drop=FALSE]%*%y),truth,tolerance=1e-10)
        for(rhs in list(y,cbind(first=y,second=2*y))) {
          a <- p07_capture(npreghat(bw,txdat=x,exdat=ex,s=1,y=rhs,output="apply"))
          expect_length(a$warnings,1L)
          expect_equal(as.vector(a$value),as.vector(H%*%rhs),tolerance=1e-10)
          expect_null(attr(a$value,".np.empty.rows",exact=TRUE))
        }
        z <- p07_capture(npreghat(bw,txdat=x,exdat=ex,s=1,y=y,output="constraint"))
        expect_length(z$warnings,1L)
        expect_equal(as.vector(z$value),as.vector(t(H)*y),tolerance=1e-10)
        pred <- p07_capture(predict(H,newdata=ex,y=y,output="apply"))
        expect_length(pred$warnings,1L)
        expect_equal(as.vector(pred$value),as.vector(H%*%y),tolerance=1e-10)
        expect_error(npreghat(bw,txdat=x,exdat=ex,s=1,.np.require.finite=TRUE))
      }
})

test_that("first hats preserve compact empty rows and strict internal helpers", {
  old <- options(np.messages=FALSE,np.tree=FALSE,matprod=getOption("matprod"),
    np.npreghat.apply.memory.threshold.mb=getOption("np.npreghat.apply.memory.threshold.mb"))
  on.exit(options(old),add=TRUE)
  ns <- asNamespace(environmentName(environment(npregbw)))
  x <- data.frame(x=seq(-1,1,length.out=41)); ex <- data.frame(x=c(.1,9,.3))
  y <- x$x^2
  bw <- npregbw(xdat=x,ydat=y,bws=.6,regtype="lp",degree=2,
                 ckertype="epanechnikov",bandwidth.compute=FALSE)
  for(tree in list(FALSE,TRUE,"auto")) for(mp in c("default","internal")) {
    options(np.tree=tree,matprod=mp)
    h <- p07_capture(npreghat(bw,txdat=x,exdat=ex,s=1))
    expect_length(h$warnings,1L)
    expect_true(all(is.na(h$value[2,])))
    expect_equal(drop(h$value[c(1,3),]%*%y),2*ex$x[c(1,3)],tolerance=1e-10)
    expect_error(get(".npreghat_complete",ns)(bw,txdat=x,exdat=ex,s=1))
  }
  expect_error(get(".npreghat_exact_lp_matrix_from_kernel_weights",ns)(
    bw,x,ex,s=1,degree=2))
  expect_error(npreghat(bw,txdat=x,exdat=ex,s=1,ridge=.1),"nonzero 'ridge'")
  expect_error(get(".npreghat_solve_eval",ns)(
    cbind(1,x$x,x$x^2),c(0,1,.2),rep(Inf,nrow(x)),0,
    canonical.lp=TRUE,allow.empty.rows=TRUE),"non-finite system")
  # Beta has a separate native partial-row policy: preserve its supported
  # direct multiresponse path, without admitting partial results there.
  xx <- data.frame(x=seq(.1,.9,length.out=41),
                   f=factor(rep(c("a","b"),length.out=41),levels=c("a","b","c")))
  ee <- data.frame(x=c(.25,.25),f=factor(c("a","c"),levels=levels(xx$f)))
  bb <- npregbw(xdat=xx,ydat=xx$x^2,bws=c(.18,0),regtype="lp",degree=2,
    ckertype="beta",ckerbound="fixed",ckerlb=0,ckerub=1,bandwidth.compute=FALSE)
  options(np.npreghat.apply.memory.threshold.mb=0)
  out <- p07_capture(npreghat(bb,txdat=xx,exdat=ee[1,,drop=FALSE],s=1,
                    y=cbind(xx$x^2,2*xx$x^2),output="apply"))
  expect_length(out$warnings,0L)
  expect_equal(as.vector(out$value[1,]),c(.5,1),tolerance=1e-10)
  expect_error(npreghat(bb,txdat=xx,exdat=ee,s=1,
                    y=cbind(xx$x^2,2*xx$x^2),output="apply"))
})
