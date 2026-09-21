.beta_quad_factor <- function(query, donor, h, order=2L, lower=0, upper=1,
                              integration=FALSE) {
  width <- upper-lower
  coef <- (-1)^(seq_len(order/2L)+1L)*choose(order/2L,seq_len(order/2L))
  ans <- numeric(length(donor))
  for(s in seq_along(coef)) {
    tau <- (width/h)^2/s
    value <- dbeta((donor-lower)/width,1+(query-lower)/width*tau,
                   1+(upper-query)/width*tau)/width
    if(integration)value[(donor==lower | donor==upper) & tau>0] <- 0
    ans <- ans+coef[s]*value
  }
  ans
}
.beta_quad_radii <- function(x,q,k,type) {
  if(type=="fixed")return(rep(k,length(x)))
  if(type=="generalized_nn")return(rep(sort(abs(x-q))[k],length(x)))
  vapply(seq_along(x),function(j)sort(abs(x[-j]-x[j]))[k],numeric(1))
}
test_that("beta PDF CVLS integrates boundary donors by their interior limit", {
  old <- options(np.messages=FALSE,np.largeh=FALSE,np.largelambda=FALSE,np.tree=FALSE)
  on.exit(options(old),add=TRUE)
  unit <- c(0,0,.13,.23,.38,.51,.67,.8,1,1); n <- length(unit)
  for(width in c(1,7)) for(order in c(2L,4L,6L,8L))
    for(type in c("fixed","generalized_nn","adaptive_nn")) {
      lo <- -2; hi <- lo+width; x <- lo+width*unit
      h <- if(type=="fixed").3*width else 4L
      b <- npudensbw(dat=data.frame(x),bws=h,bwmethod="cv.ls",
        bwtype=type,ckertype="beta",ckerorder=order,
        ckerbound="fixed",ckerlb=lo,ckerub=hi,bandwidth.compute=FALSE)
      query <- seq(lo,hi,length.out=81)
      weights <- rep(width/80,81);weights[c(1,81)] <- weights[c(1,81)]/2
      fit <- vapply(query,function(q)mean(.beta_quad_factor(q,x,
        .beta_quad_radii(x,q,h,type),order,lo,hi,TRUE)),numeric(1))
      cross <- vapply(seq_along(x),function(i)mean(.beta_quad_factor(x[i],x[-i],
        .beta_quad_radii(x[-i],x[i],h,type),order,lo,hi)),numeric(1))
      actual <- getFromNamespace("npudensbw.bandwidth","np")(dat=data.frame(x),
        bws=b,bandwidth.compute=TRUE,eval.only=TRUE,nmulti=1L)$fval
      expect_equal(as.numeric(actual),2*mean(cross)-sum(weights*fit^2),
        tolerance=2e-10,info=paste(width,order,type))
      point <- fitted(npudens(bws=b,tdat=data.frame(x),edat=data.frame(x=c(lo,hi))))
      expected <- vapply(c(lo,hi),function(q)mean(.beta_quad_factor(q,x,
        .beta_quad_radii(x,q,h,type),order,lo,hi)),numeric(1))
      expect_equal(as.numeric(point),expected,tolerance=2e-11)
    }
})
test_that("zero-concentration beta components retain their uniform integral", {
  old <- options(np.messages=FALSE,np.largeh=FALSE,np.largelambda=FALSE)
  on.exit(options(old),add=TRUE)
  dat <- data.frame(x=c(0,0,.2,.4,.7,1,1))
  for(order in c(2L,4L,6L,8L)) {
    b <- npudensbw(dat=dat,bws=1e200,bwmethod="cv.ls",ckertype="beta",
      ckerorder=order,ckerbound="fixed",ckerlb=0,ckerub=1,bandwidth.compute=FALSE)
    actual <- getFromNamespace("npudensbw.bandwidth","np")(dat=dat,bws=b,
      bandwidth.compute=TRUE,eval.only=TRUE,nmulti=1L)$fval
    expect_equal(as.numeric(actual),1,tolerance=2e-12)
    expect_equal(as.numeric(fitted(npudens(bws=b,tdat=dat))),rep(1,nrow(dat)),
      tolerance=2e-12)
  }
})
test_that("conditional beta PDF quadrature shares the representative in every fold", {
  old <- options(np.messages=FALSE,np.largeh=FALSE,np.largelambda=FALSE,np.tree=FALSE)
  on.exit(options(old),add=TRUE)
  x <- c(.03,.09,.17,.23,.36,.41,.58,.69,.74,.83,.91,.98)
  ys <- list(c(0,0,.13,.23,.38,.51,.67,.8,.85,.92,1,1),
             c(.13,1,.51,0,.85,.23,1,.67,.38,0,.8,.92))
  n <- length(x)
  for(p in c(1L,2L)) for(type in c("fixed","generalized_nn","adaptive_nn"))
    for(reg in c("lc","ll","lp")) {
      y <- as.data.frame(ys[seq_len(p)]); names(y) <- paste0("y",seq_len(p))
      h <- if(type=="fixed").35 else 5L
      args <- list(xdat=data.frame(x),ydat=y,bws=rep(h,p+1L),
        bwmethod="cv.ls",bwtype=type,regtype=reg,cykertype="beta",
        cykerbound="fixed",cykerlb=rep(0,p),cykerub=rep(1,p),
        bandwidth.compute=FALSE,cvls.quadrature.grid="uniform",
        cvls.quadrature.points=c(81L,31L))
      degree <- switch(reg,lc=0L,ll=1L,lp=2L)
      if(reg=="lp")args$degree <- degree
      b <- do.call(npcdensbw,args)
      q <- if(p==1L)81L else 31L
      g <- seq(0,1,length.out=q);w <- rep(1/(q-1),q);w[c(1,q)] <- w[c(1,q)]/2
      ids <- as.matrix(expand.grid(rep(list(seq_len(q)),p)))
      grid <- matrix(g[ids],ncol=p); quadw <- apply(matrix(w[ids],ncol=p),1,prod)
      loss <- vapply(seq_len(n),function(i) {
        donors <- setdiff(seq_len(n),i)
        hx <- .beta_quad_radii(x[donors],x[i],h,type)
        wx <- dnorm((x[i]-x[donors])/hx)/hx
        design <- outer(x[donors]-x[i],0:degree,"^")
        influence <- solve(crossprod(design,wx*design),t(design*wx))[1,]
        yy <- y[donors,,drop=FALSE]
        one <- function(query,integration) {
          factors <- rep(1,length(donors))
          for(d in seq_len(p))factors <- factors*.beta_quad_factor(query[d],yy[[d]],
            .beta_quad_radii(yy[[d]],query[d],h,type),integration=integration)
          sum(influence*factors)
        }
        full <- apply(grid,1,one,integration=TRUE)
        sum(quadw*full^2)-2*one(as.numeric(y[i,]),FALSE)
      },numeric(1))
      actual <- getFromNamespace(".npcdensbw_eval_only","np")(data.frame(x),y,b)$objective
      expect_equal(as.numeric(actual),-mean(loss),tolerance=3e-10,
        info=paste(p,type,reg))
    }
})
