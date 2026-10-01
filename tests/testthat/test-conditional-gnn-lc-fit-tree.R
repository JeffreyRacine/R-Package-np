gnn_lc_conditional_literal <- function(x,y,e,ey,k,kernel,cdf,training) {
  K <- function(u)if(kernel=="uniform").5*(abs(u)<1)else
    3/(4*sqrt(5))*pmax(0,1-u*u/5)
  F <- function(u) {
    z <- pmax(-1,pmin(1,u/if(kernel=="uniform")1 else sqrt(5)))
    if(kernel=="uniform")(z+1)/2 else .5+.75*z-.25*z^3
  }
  vapply(seq_len(nrow(e)),function(j) {
    w <- rep(1,nrow(x))
    for(l in seq_len(ncol(x))) {
      d <- abs(x[,l]-e[j,l]);h <- sort(if(training)d[-j]else d)[k]
      w <- w*K((x[,l]-e[j,l])/h)
    }
    d <- abs(y-ey[j]);hy <- sort(if(training)d[-j]else d)[k]
    ky <- if(cdf)F((ey[j]-y)/hy)else K((ey[j]-y)/hy)/hy
    sum(w*ky)/sum(w)
  },0.0)
}

test_that("conditional GNN fitting admission leaves search and other routes intact", {
  old <- options(np.messages=FALSE,np.tree=FALSE);on.exit(options(old))
  ns <- asNamespace("np");yes <- get("DO_TREE_YES",ns);no <- get("DO_TREE_NO",ns)
  x <- data.frame(x=seq(-1,1,length.out=40));y <- sin(seq_len(40))
  for(cdf in c(FALSE,TRUE)) {
    ctor <- if(cdf)npcdistbw else npcdensbw
    choose <- get(if(cdf)".npcdistbw_tree_code"else".npcdensbw_tree_code",ns)
    b <- ctor(xdat=x,ydat=y,bws=c(20,20),bwtype="generalized_nn",
      bwmethod="cv.ls",regtype="lc",cxkertype="epanechnikov",
      cykertype="epanechnikov",bandwidth.compute=FALSE)
    for(mode in list(FALSE,TRUE,"auto")) for(fit in c(FALSE,TRUE))
      for(engine in c("lc","lp")) for(kernel in c("epanechnikov","uniform","gaussian"))
        for(order in c(2L,4L)) for(bound in c("none","fixed")) {
          z <- b;z$regtype.engine <- engine
          z$cxkertype <- z$cykertype <- kernel;z$cxkerorder <- z$cykerorder <- order
          z$cxkerbound <- z$cykerbound <- bound;options(np.tree=mode)
          # Search and fitting have distinct, independently stated contracts.
          selected <- isTRUE(mode) || (identical(mode,"auto") && kernel!="gaussian")
          fit.allowed <- isTRUE(mode) && fit && kernel!="gaussian" &&
            order==2L && bound=="none"
          search.allowed <- !cdf && !fit && selected && bound=="none" &&
            (kernel=="epanechnikov" || (kernel%in%c("gaussian","uniform") && order==2L))
          expected <- if(fit.allowed || search.allowed)yes else no
          expect_identical(choose(z,2L,0L,fit.context=fit),expected)
        }
    options(np.tree=TRUE)
    if(cdf)expect_identical(choose(b,2L,0L,cv.context=TRUE),yes)
    for(type in c("fixed","adaptive_nn")) {
      z <- b;z$type <- type
      expect_identical(choose(z,2L,0L,fit.context=FALSE),choose(z,2L,0L,fit.context=TRUE))
    }
    if(!cdf) {
      z <- b;z$method <- "cv.ml"
      expect_identical(choose(z,2L,0L,fit.context=FALSE),choose(z,2L,0L,fit.context=TRUE))
    }
  }
})

test_that("conditional GNN LC trees preserve literal kernels and uncertainty", {
  old <- options(np.messages=FALSE,np.tree=FALSE);on.exit(options(old))
  set.seed(62);n <- 192L
  x0 <- data.frame(x1=runif(n,-1,1),x2=runif(n,-1,1));x0[2,] <- x0[1,]
  y <- sin(x0$x1)+rnorm(n,sd=.5)
  for(cdf in c(FALSE,TRUE))for(p in 1:2)for(kernel in c("epanechnikov","uniform"))
    for(external in c(FALSE,TRUE)) {
      x <- x0[,seq_len(p),drop=FALSE];k <- 125L
      e <- if(external).99*x[seq(3,n,by=11),,drop=FALSE]else x
      ey <- if(external)y[seq(3,n,by=11)]+.001 else y
      ctor <- if(cdf)npcdistbw else npcdensbw;fit <- if(cdf)npcdist else npcdens
      b <- ctor(xdat=x,ydat=y,bws=rep(k,p+1),bwtype="generalized_nn",bwmethod="cv.ls",
        regtype="lc",cxkertype=kernel,cykertype=kernel,bandwidth.compute=FALSE)
      args <- list(bws=b,txdat=x,tydat=y,gradients=TRUE,se=TRUE)
      if(external){args$exdat <- e;args$eydat <- ey}
      options(np.tree=FALSE);off <- do.call(fit,args)
      options(np.tree=TRUE);on <- do.call(fit,args)
      ref <- gnn_lc_conditional_literal(x,y,e,ey,k,kernel,cdf,!external)
      expect_true(max(abs(fitted(off)-ref))<2e-9)
      expect_true(max(abs(fitted(on)-ref))<2e-9)
      expect_true(max(abs(fitted(on)-fitted(off)))<2e-9)
      expect_true(max(abs(se(on)-se(off)))<2e-9)
      expect_true(max(abs(gradients(on)-gradients(off)))<2e-9)
      expect_true(max(abs(gradients(on,se=TRUE)-gradients(off,se=TRUE)))<2e-9)
    }
})

test_that("mixed conditional GNN LC fits and formula predictions keep support and order", {
  old <- options(np.messages=FALSE,np.tree=FALSE);on.exit(options(old))
  set.seed(47);n <- 160L
  d <- data.frame(x=runif(n,-1,1),u=factor(rep(c("a","b"),length.out=n)),
    o=ordered(rep(c(.5,1,2,4),length.out=n)))
  d$y <- sin(d$x)+rnorm(n,sd=.5)
  e <- d[seq(3,n,by=9),]
  for(cdf in c(FALSE,TRUE)) {
    ctor <- if(cdf)npcdistbw else npcdensbw;fit <- if(cdf)npcdist else npcdens
    b <- ctor(y~x+u+o,data=d,bws=c(110,110,.2,.3),bwtype="generalized_nn",
      bwmethod="cv.ls",regtype="lc",cxkertype="epanechnikov",cykertype="epanechnikov",
      oxkertype="racineliyan",bandwidth.compute=FALSE)
    output <- lapply(c(FALSE,TRUE,FALSE,TRUE),function(tree) {
      options(np.tree=tree);g <- fit(bws=b,gradients=TRUE,se=TRUE)
      c(fitted(g),se(g),gradients(g),gradients(g,se=TRUE),predict(g,newdata=e))
    })
    for(j in 2:4) {
      expect_identical(is.na(output[[1]]),is.na(output[[j]]))
      expect_true(max(abs(output[[1]]-output[[j]]),na.rm=TRUE)<2e-9)
    }
  }
})
