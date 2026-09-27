# Independent GNN fit oracle: raw polynomial WLS and literal NN distances.
gnn_fit_literal <- function(x, y, e, k, degree, basis, kernel, training) {
  x <- as.matrix(x); e <- as.matrix(e); p <- ncol(x); n <- nrow(x)
  terms <- as.matrix(expand.grid(rep(list(0:degree), p)))
  if (basis == "glp") terms <- terms[rowSums(terms) <= degree,,drop=FALSE]
  if (basis == "additive") terms <- terms[rowSums(terms > 0) <= 1,,drop=FALSE]
  design <- function(z) vapply(seq_len(nrow(terms)), function(a)
    apply(sweep(z, 2, terms[a,], "^"), 1, prod), numeric(nrow(z)))
  B <- design(x); E <- design(e)
  K <- if(kernel == "uniform") function(u) .5*(abs(u)<1) else
    function(u) 3/(4*sqrt(5))*pmax(0,1-u*u/5)
  vapply(seq_len(nrow(e)), function(j) {
    w <- rep(1,n)
    for(l in seq_len(p)) {
      d <- abs(x[,l]-e[j,l]); h <- sort(if(training)d[-j]else d)[k[l]]
      w <- w*K((x[,l]-e[j,l])/h)
    }
    sum(E[j,]*solve(crossprod(B,w*B),crossprod(B,w*y)))
  }, 0.0)
}

test_that("GNN fitting admission is explicit, compact and route limited", {
  old <- options(np.messages=FALSE,np.tree=FALSE); on.exit(options(old))
  ns <- asNamespace("npRmpi"); choose <- get(".npreg_fit_tree_code",ns)
  b <- npregbw(xdat=data.frame(x=seq(-1,1,length.out=40)),ydat=sin(1:40),
    bws=20,bwtype="generalized_nn",regtype="lp",degree=2L,
    ckertype="epanechnikov",bandwidth.compute=FALSE)
  yes <- get("DO_TREE_YES",ns); no <- get("DO_TREE_NO",ns)
  for(mode in list(FALSE,TRUE,"auto")) for(type in c("fixed","generalized_nn","adaptive_nn"))
    for(kernel in c("epanechnikov","uniform","gaussian")) for(order in c(2L,4L))
      for(bound in c("none","fixed")) {
        z <- b; z$type <- type; z$ckertype <- kernel
        z$ckerorder <- order; z$ckerbound <- bound; options(np.tree=mode)
        expected <- if(isTRUE(mode) && type %in% c("generalized_nn","adaptive_nn") &&
          kernel %in% c("epanechnikov","uniform") && order==2L && bound=="none")yes else no
        expect_identical(choose(z,1L,0L),expected,
          info=paste(mode,type,kernel,order,bound))
      }
  # LC retains the generic selector; LP0 retains its LP engine identity.
  for(regtype in c("lc","lp")) {
    z <- npregbw(xdat=data.frame(x=seq_len(40)),ydat=sin(1:40),
      bws=20,bwtype="generalized_nn",regtype=regtype,degree=0L,
      ckertype="epanechnikov",bandwidth.compute=FALSE)
    for(mode in list(FALSE,TRUE,"auto")) {
      options(np.tree=mode)
      expect_identical(choose(z,1L,0L),
        if(regtype == "lc" || isTRUE(mode))
          get("npDoTreeOrCategoricalCompress",ns)(ncon=1L,ncat=0L,bws=z) else no)
    }
  }
})

test_that("GNN tree fitting preserves literal WLS and uncertainty contracts", {
  old <- options(np.messages=FALSE,np.tree=FALSE); on.exit(options(old))
  set.seed(20260925); n <- 192L
  dat <- as.data.frame(matrix(runif(n*4,-1,1),n,4)); names(dat) <- paste0("x",1:4)
  dat[2,] <- dat[1,]; y <- sin(2*dat$x1)+.4*dat$x2^2+rnorm(n,sd=.12)
  cases <- data.frame(p=c(1,1,2,2,3,4),d=c(1,3,1,2,2,1),
    basis=c("glp","glp","glp","tensor","additive","glp"))
  for(a in seq_len(nrow(cases))) for(kernel in c("epanechnikov","uniform"))
    for(bernstein in c(FALSE,TRUE)) for(external in c(FALSE,TRUE)) {
      q <- cases[a,]; x <- dat[,seq_len(q$p),drop=FALSE]; k <- rep(145L,q$p)
      e <- if(external)x[seq(5,n,length.out=23),,drop=FALSE]+.001 else x
      options(np.tree=FALSE)
      b <- npregbw(xdat=x,ydat=y,bws=k,bwtype="generalized_nn",regtype="lp",
        degree=rep(q$d,q$p),basis=q$basis,bernstein.basis=bernstein,
        ckertype=kernel,bandwidth.compute=FALSE)
      args <- list(bws=b,txdat=x,tydat=y,gradients=TRUE,se=TRUE)
      if(external)args$exdat <- e
      off <- do.call(npreg,args); options(np.tree=TRUE); on <- do.call(npreg,args)
      ref <- gnn_fit_literal(x,y,e,k,q$d,q$basis,kernel,!external)
      check <- function(v,w,what) expect_true(max(abs(v-w))<2e-9,
        info=paste(a,kernel,bernstein,external,what,max(abs(v-w))))
      check(fitted(off),ref,"off literal");check(fitted(on),ref,"on literal")
      check(fitted(on),fitted(off),"point");check(se(on),se(off),"SE")
      check(gradients(on),gradients(off),"gradient")
      check(gradients(on,se=TRUE),gradients(off,se=TRUE),"gradient SE")
    }
})

test_that("GNN fitting consumers preserve mixed-data operators and predictions", {
  old <- options(np.messages=FALSE,np.tree=FALSE); on.exit(options(old))
  ns <- asNamespace("npRmpi")
  set.seed(813); n <- 144L
  x <- data.frame(x1=runif(n,-1,1), x2=runif(n,-1,1),
    u=factor(rep(c("a","b","c"),length.out=n)),o=ordered(rep(1:3,length.out=n)))
  y <- sin(x$x1)+x$x2^2+as.integer(x$u)/5+rnorm(n,sd=.1)
  e <- x[seq(3,n,by=9),]; e$x1 <- .99*e$x1
  for(degree in 1:2) {
    b <- npregbw(xdat=x,ydat=y,bws=c(115,115,.2,.3),bwtype="generalized_nn",
      regtype="lp",degree=rep(degree,2),bernstein.basis=TRUE,
      ckertype="epanechnikov",okertype="liracine",bandwidth.compute=FALSE)
    result <- lapply(c(FALSE,TRUE),function(tree) {
      options(np.tree=tree)
      fit <- npreg(bws=b,txdat=x,tydat=y,exdat=e,gradients=TRUE,se=TRUE)
      direct <- get(".np_regression_direct",ns)(bws=b,txdat=x,tydat=y,exdat=e)
      H <- as.matrix(npreghat(bws=b,txdat=x,exdat=e,output="matrix"))
      applied <- as.numeric(npreghat(bws=b,txdat=x,exdat=e,y=y,output="apply"))
      loo <- as.numeric(npreghat(bws=b,txdat=x,y=y,leave.one.out=TRUE,output="apply"))
      expect_true(max(abs(drop(H%*%y)-applied))<2e-9)
      c(fitted(fit),se(fit),gradients(fit),gradients(fit,se=TRUE),
        unlist(direct),H,applied,loo)
    })
    expect_true(max(abs(result[[1]]-result[[2]]))<2e-9)
  }
  d <- cbind(y=y,x)
  options(np.tree=FALSE)
  b <- npregbw(y~x1+x2+u+o,data=d,bws=c(115,115,.2,.3),
    bwtype="generalized_nn",regtype="lp",degree=c(2L,2L),
    bernstein.basis=TRUE,ckertype="epanechnikov",bandwidth.compute=FALSE)
  g <- npreg(bws=b); off <- predict(g,newdata=e)
  options(np.tree=TRUE); on <- predict(g,newdata=e)
  expect_true(max(abs(on-off))<2e-9)
})

test_that("GNN paired conditional uncertainty retains its X-weight contract", {
  old <- options(np.messages=FALSE,np.tree=FALSE); on.exit(options(old))
  set.seed(41); n <- 96L
  x <- data.frame(x=runif(n,-1,1)); y <- data.frame(y=sin(x$x)+rnorm(n,sd=.3))
  b <- npcdensbw(xdat=x,ydat=y,bws=c(65,65),bwtype="generalized_nn",
    regtype="lp",degree=2L,bernstein.basis=TRUE,
    cxkertype="epanechnikov",cykertype="epanechnikov",bandwidth.compute=FALSE)
  pair <- getFromNamespace(".np_conditional_lp_pair_se","npRmpi")
  upper <- x[seq(4,n,by=8),,drop=FALSE]; lower <- upper-.01
  ey <- y[seq(4,n,by=8),,drop=FALSE]
  for(cdf in c(FALSE,TRUE)) {
    options(np.tree=FALSE); off <- pair(b,x,y,upper,lower,ey,cdf=cdf)
    options(np.tree=TRUE); on <- pair(b,x,y,upper,lower,ey,cdf=cdf)
    expect_identical(is.na(on),is.na(off))
    expect_true(max(abs(on-off),na.rm=TRUE)<2e-9)
  }
})


test_that("GNN full-support fitting preserves weights and uncertainty", {
  old <- options(np.messages=FALSE,np.tree=FALSE,np.extendednn=TRUE)
  on.exit(options(old))
  n <- 64L; x <- data.frame(x=seq(-1,1,length.out=n))
  y <- sin(3*x$x)+x$x^2
  for(kernel in c("epanechnikov","uniform"))for(k in c(n-2L,n-1L,n+2L))
    for(external in c(FALSE,TRUE)) {
      b <- npregbw(xdat=x,ydat=y,bws=k,bwtype="generalized_nn",
        regtype="lp",degree=2L,bernstein.basis=TRUE,
        ckertype=kernel,bandwidth.compute=FALSE)
      args <- list(bws=b,txdat=x,tydat=y,se=TRUE,gradients=TRUE)
      if(external)args$exdat <- data.frame(x=c(-1,0,1))
      z <- lapply(c(FALSE,TRUE),function(tree) {
        options(np.tree=tree); g <- do.call(npreg,args)
        c(fitted(g),se(g),gradients(g),gradients(g,se=TRUE))
      })
      expect_true(max(abs(z[[1]]-z[[2]]))<2e-9)
    }
})

test_that("GNN tree apply preserves batched RHS donor identities", {
  old <- options(np.messages=FALSE,np.tree=FALSE); on.exit(options(old))
  set.seed(813); n <- 144L
  x <- data.frame(x1=runif(n,-1,1),x2=runif(n,-1,1))
  Y <- cbind(sin(x$x1)+x$x2^2,cos(x$x1)-x$x2,x$x1*x$x2)
  for(permute in c(FALSE,TRUE)) {
    ix <- if(permute)sample.int(n) else seq_len(n)
    xx <- x[ix,]; yy <- Y[ix,,drop=FALSE]; saved <- yy
    b <- npregbw(xdat=xx,ydat=yy[,1],bws=c(115,115),
      bwtype="generalized_nn",regtype="lp",degree=c(2L,2L),
      bernstein.basis=TRUE,ckertype="epanechnikov",bandwidth.compute=FALSE)
    for(loo in c(FALSE,TRUE)) for(external in c(FALSE,TRUE)) {
      if(loo && external)next
      args <- list(bws=b,txdat=xx,leave.one.out=loo)
      if(external)args$exdat <- x[seq(3,n,by=9),]
      outputs <- lapply(c(FALSE,TRUE),function(tree){
        options(np.tree=tree)
        H <- do.call(npreghat,c(args,list(output="matrix")))
        A <- do.call(npreghat,c(args,list(y=yy,output="apply")))
        expect_true(max(abs(A-H%*%yy))<2e-9)
        expect_identical(yy,saved)
        A
      })
      expect_true(max(abs(outputs[[1]]-outputs[[2]]))<2e-9)
    }
  }
})
