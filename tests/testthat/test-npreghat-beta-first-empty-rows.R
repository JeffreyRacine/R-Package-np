p08_capture <- function(expr) {
  observed <- new.env(parent=emptyenv())
  observed$warnings <- character()
  value <- withCallingHandlers(expr, warning=function(w) {
    observed$warnings <- c(observed$warnings,conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  list(value=value,warnings=observed$warnings)
}

test_that("beta first LP hats preserve supported rows and mark empty rows", {
  skip_on_cran()
  old <- options(np.messages=FALSE,np.tree=FALSE)
  on.exit(options(old),add=TRUE)
  x <- data.frame(x=seq(.1,.9,length.out=33),
    f=factor(rep(c("a","b"),length.out=33),levels=c("a","b","c")))
  ex <- data.frame(x=c(.25,.25,.65),f=factor(c("a","c","b"),levels=levels(x$f)))
  for(type in c("fixed","generalized_nn","adaptive_nn"))
    for(spec in list(list(regtype="ll"),list(regtype="lp",degree=1),
                    list(regtype="lp",degree=2),
                    list(regtype="lp",degree=2,bernstein.basis=TRUE))) {
      y <- if(identical(spec$degree,2))x$x^2 else 2+3*x$x
      b <- do.call(npregbw,c(list(xdat=x,ydat=y,
        bws=c(if(type=="fixed").18 else 23,0),bwtype=type,
        ckertype="beta",ckerbound="fixed",ckerlb=0,ckerub=1,
        bandwidth.compute=FALSE),spec))
      h <- p08_capture(npreghat(b,txdat=x,exdat=ex,s=1))
      H <- h$value
      expect_length(h$warnings,1L)
      expect_match(h$warnings,"all computed kernel weights are zero")
      expect_true(all(is.na(H[2,])))
      expect_true(all(is.finite(H[c(1,3),])))
      expect_true(is.na(attr(H,"ridge.used")[2]))
      expect_null(attr(H,".np.empty.rows",exact=TRUE))
      truth <- if(identical(spec$degree,2))2*ex$x[c(1,3)] else c(3,3)
      expect_equal(drop(H[c(1,3),,drop=FALSE]%*%y),truth,tolerance=1e-10)
      for(rhs in list(y,cbind(y,2*y))) {
        a <- p08_capture(npreghat(b,txdat=x,exdat=ex,s=1,y=rhs,output="apply"))
        expect_length(a$warnings,1L)
        expect_equal(as.vector(a$value),as.vector(H%*%rhs),tolerance=1e-10)
        expect_null(attr(a$value,".np.empty.rows",exact=TRUE))
      }
      z <- p08_capture(npreghat(b,txdat=x,exdat=ex,s=1,y=y,output="constraint"))
      expect_length(z$warnings,1L)
      expect_equal(as.vector(z$value),as.vector(t(H)*y),tolerance=1e-10)
      p <- p08_capture(predict(H,newdata=ex,y=y,output="apply"))
      expect_length(p$warnings,1L)
      expect_equal(as.vector(p$value),as.vector(H%*%y),tolerance=1e-10)
      expect_error(npreghat(b,txdat=x,exdat=ex,s=1,.np.require.finite=TRUE))
    }
})

test_that("beta blocked apply keeps empty row indices across block boundaries", {
  skip_on_cran()
  old <- options(np.messages=FALSE,np.tree=FALSE,
    np.npreghat.apply.memory.threshold.mb=0,matprod=getOption("matprod"))
  on.exit(options(old),add=TRUE)
  x <- data.frame(x=seq(.1,.9,length.out=33),
    f=factor(rep(c("a","b"),length.out=33),levels=c("a","b","c")))
  ex <- data.frame(x=seq(.2,.8,length.out=259),
    f=factor(rep("a",259),levels=levels(x$f)))
  bad <- c(1L,128L,256L,257L,259L);ex$f[bad] <- "c"
  b <- npregbw(xdat=x,ydat=x$x^2,bws=c(.18,0),regtype="lp",degree=2,
    ckertype="beta",ckerbound="fixed",ckerlb=0,ckerub=1,bandwidth.compute=FALSE)
  rhs <- cbind(x$x^2,2*x$x^2,1+x$x)
  for(mp in c("default","internal")) {
    options(matprod=mp)
    a <- p08_capture(npreghat(b,txdat=x,exdat=ex,s=1,y=rhs,output="apply"))
    expect_length(a$warnings,1L)
    expect_true(all(is.na(a$value[bad,])))
    expect_equal(unname(a$value[-bad,]),
      cbind(2*ex$x[-bad],4*ex$x[-bad],1),tolerance=1e-10)
  }
  expect_error(npreghat(b,txdat=x,exdat=ex,s=1,y=rhs,output="apply",
    .np.require.finite=TRUE))
  # Other native partial-row contracts are not widened by this repair.
  expect_error(npreghat(b,txdat=x,exdat=ex,s=0))
  expect_error(npreghat(b,txdat=x,exdat=ex,s=1,ridge=.1),"nonzero 'ridge'")
})
