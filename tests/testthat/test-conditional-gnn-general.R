cgnn_general_evaluate <- function(x,y,b) get('.npcdensbw_eval_only',asNamespace('np'))(x,y,b,invalid.penalty='dbmax')$objective

# Independent q-space partition, then exact polynomial quadrature. Within a
# NN/support interval h(q)=scale*abs(q-anchor). With u=1/(q-anchor), the squared
# fit times dq is a polynomial of degree <=16 for the compact kernels. Nine
# Gauss-Legendre nodes integrate it exactly in real arithmetic. Eleven nodes
# provide a separate floating-arithmetic discrepancy check, not a new tolerance.
cgnn_compact_integral <- function(values,k,family,order,weights=NULL,pair=NULL) {
  rule <- function(n) {
    b <- seq_len(n-1)/sqrt(4*seq_len(n-1)^2-1)
    J <- matrix(0,n,n);J[cbind(1:(n-1),2:n)] <- b
    J <- J+t(J);e <- eigen(J,symmetric=TRUE)
    list(x=e$values,w=2*e$vectors[1,]^2)
  }
  rules <- list(rule(9L),rule(11L))
  cuts <- cgnn_general_cuts(values,family,order,k)
  index <- min(k,length(values)-1L);scale <- k/index
  total <- c(0,0);magnitude <- 0
  for(at in seq_len(length(cuts)-1L)) {
    lo <- cuts[at];hi <- cuts[at+1L]
    mid <- if(!is.finite(lo))hi-max(1,abs(hi))else
      if(!is.finite(hi))lo+max(1,abs(lo))else lo/2+hi/2
    anchor <- values[order(abs(mid-values))[index]]
    ends <- sort(c(1/(lo-anchor),1/(hi-anchor)))
    half <- (ends[2]-ends[1])/2;center <- ends[1]/2+ends[2]/2
    stopifnot(is.finite(half),is.finite(center),half>=0)
    for(r in 1:2) {
      u <- center+half*rules[[r]]$x
      val <- vapply(u,function(v) {
        z <- cgnn_general_kernel((1+(anchor-values)*v)/scale,family,order)
        if(is.null(pair))sum(weights*z)^2/scale^2 else prod(z[pair])/scale^2
      },0)
      total[r] <- total[r]+half*sum(rules[[r]]$w*val)
      magnitude <- magnitude+half*sum(abs(rules[[r]]$w*val))
    }
  }
  c(value=total[1],error=abs(total[1]-total[2])+64*.Machine$double.eps*magnitude)
}

# Independent literal deleted WLS / coordinate-integral oracle. No package
# radius, overlap, basis, or influence helper is used by this reference.
cgnn_general_kernel <- function(z, family, order=2L) {
  v <- z*z
  if(family == 'uniform') return(.5*(abs(z)<1))
  if(family == 'gaussian') return(dnorm(z)*switch(as.character(order),
    '2'=1, '4'=1.5-.5*v, '6'=1.875+v*(-1.25+.125*v),
    '8'=2.1875+v*(-2.1875+v*(.4375-.02083333333*v))))
  value <- switch(as.character(order),
    '2'=.33541019662496845446*(1-v/5),
    '4'=.008385254916*(-15+7*v)*(-5+v),
    '6'=.33541019662496845446*(2.734375+v*(-3.28125+.721875*v))*(1-.2*v),
    '8'=.33541019662496845446*(3.5888671875+v*(-7.8955078125+
       v*(4.1056640625-.5865234375*v)))*(1-.2*v))
  value[v>=5] <- 0
  value
}

cgnn_general_radius <- function(z, values, k) {
  index <- min(k,length(values)-1L)
  sort(abs(values-z),partial=index)[index]*(k/index)
}

cgnn_general_cuts <- function(y, family, order, k) {
  cuts <- c(-Inf, y, as.vector(outer(y,y,'+')/2), Inf)
  if(family != 'gaussian') {
    radius <- (if(family=='uniform')1 else sqrt(5))*(k/min(k,length(y)-1L))
    for(sign in c(-1,1)) if(1-sign*radius!=0)
      cuts <- c(cuts,as.vector(outer(y,sign*radius*y,'-')/(1-sign*radius)))
  }
  sort(unique(cuts))
}

cgnn_general_reference <- function(x, y, kx, ky, degree=0L,
                                   family='gaussian', order=2L) {
  n <- length(y); basis <- outer(x,0:degree,'^')
  result <- matrix(0,n,3)
  for(i in seq_len(n)) {
    xx <- x[-i]; yy <- y[-i]; A <- basis[-i,,drop=FALSE]
    w <- dnorm((x[i]-xx)/cgnn_general_radius(x[i],xx,kx))
    weights <- w*as.vector(A%*%solve(crossprod(A,w*A),basis[i,]))
    density <- function(q) vapply(q,function(z) {
      h <- cgnn_general_radius(z,yy,ky)
      sum(weights*cgnn_general_kernel((z-yy)/h,family,order)/h)
    },0)
    if(family != 'gaussian') {
      z <- cgnn_compact_integral(yy,ky,family,order,weights)
      result[i,] <- c(z['value'],density(y[i]),z['error'])
    } else {
      cuts <- cgnn_general_cuts(yy,family,order,ky)
      ints <- vapply(seq_len(length(cuts)-1L),function(j) {
        z <- integrate(function(q)density(q)^2,cuts[j],cuts[j+1L],
           abs.tol=1e-13/length(cuts),rel.tol=1e-11,subdivisions=1000L)
        c(z$value,z$abs.error)
      },c(0,0))
      result[i,] <- c(sum(ints[1,]),density(y[i]),sum(ints[2,]))
    }
  }
  z <- colMeans(result)
  c(I1=z[1],I2=z[2],score=2*z[2]-z[1],error=z[3])
}

test_that('general conditional GNN kernels use the literal deleted fit', {
  old <- options(np.messages=FALSE,np.tree=FALSE)
  on.exit(options(old),add=TRUE)
  x <- c(-.67932861763984,.597772665787488,-.135870280675590,-.180028952658176,
    -.0812444845214486,-.0922158942557871,-.119163474533707,.565201563760638,
    -.167390177026391,-.617473398335278,-.452014316339046,.213771468494087)
  y <- c(.117393371024781,-.215783902381389,-.752371990755028,1.55294725882662,
    .776713307792541,.329678953352629,-1.89405512568533,-1.99349412668255,
    -.220072074888336,.198501257716297,.563714301921436,-.880794421005045)
  for(family in c('gaussian','epanechnikov','uniform'))
    for(order in if(family=='uniform')2L else c(2L,4L,6L,8L))
      for(degree in c(0L,2L)) {
        ref <- cgnn_general_reference(x,y,8L,3L,degree,family,order)
        expect_lt(ref['error'],2.5e-10)
        ctl <- if(degree==0L)list(regtype='lc')else
          list(regtype='lp',degree=2L,bernstein.basis=TRUE)
        b <- do.call(npcdensbw,c(list(xdat=data.frame(x=x),ydat=y,bws=c(3,8),
          bwtype='generalized_nn',bwmethod='cv.ls',cxkertype='gaussian',
          cykertype=family,bandwidth.compute=FALSE),ctl,
          if(family=="uniform")list()else list(cykerorder=order)))
        observed <- cgnn_general_evaluate(data.frame(x=x),y,b)
        expect_true(abs(observed-ref['score']) <= 1e-9,
          info=paste(family,order,degree))
      }
})

# Reuse the independently written kernel/radius definitions from this file.
# The tensor oracle below integrates each deleted donor pair coordinate by
# coordinate, then explicitly sums the tiny pair matrix. Only tests store it.
cgnn_general_tensor_reference <- function(x,y,kx,ky,degree=0L,
                                          family='gaussian',category=NULL) {
  y <- as.matrix(y); n <- nrow(y); dims <- ncol(y)
  B <- outer(x,0:degree,'^'); answer <- matrix(0,n,3)
  for(i in seq_len(n)) {
    keep <- seq_len(n)!=i; yy <- y[keep,,drop=FALSE]
    A <- B[keep,,drop=FALSE]
    w <- dnorm((x[i]-x[keep])/cgnn_general_radius(x[i],x[keep],kx))
    w <- w*as.vector(A%*%solve(crossprod(A,w*A),B[i,]))
    overlap <- matrix(1,n-1,n-1); bound <- matrix(0,n-1,n-1)
    cross <- rep(1,n-1)
    for(d in seq_len(dims)) {
      values <- yy[,d]
      cuts <- cgnn_general_cuts(values,family,2L,ky)
      h <- cgnn_general_radius(y[i,d],values,ky)
      cross <- cross*cgnn_general_kernel((y[i,d]-values)/h,family)/h
      for(j in seq_len(n-1))for(k in seq_len(j)) {
        product <- function(z) vapply(z,function(q) {
          h <- cgnn_general_radius(q,values,ky)
          prod(cgnn_general_kernel((q-values[c(j,k)])/h,family)/h)
        },0)
        ints <- vapply(seq_len(length(cuts)-1L),function(at) {
          z <- integrate(product,cuts[at],cuts[at+1L],abs.tol=1e-13/length(cuts),
                         rel.tol=1e-11,subdivisions=1000L)
          c(z$value,z$abs.error)
        },c(0,0))
        val <- sum(ints[1,]); err <- sum(ints[2,]); prior <- overlap[j,k]
        bound[j,k] <- bound[j,k]*(abs(val)+err)+abs(prior)*err
        overlap[j,k] <- prior*val
        overlap[k,j] <- overlap[j,k];bound[k,j] <- bound[j,k]
      }
    }
    if(!is.null(category)) {
      # Explicit Aitchison-Aitken probability vectors over declared support.
      lev <- levels(category); lambda <- .2
      probabilities <- vapply(category[keep],function(v)
        ifelse(lev==v,1-lambda,lambda/(length(lev)-1)),numeric(length(lev)))
      cat_overlap <- crossprod(probabilities)
      overlap <- overlap*cat_overlap;bound <- bound*abs(cat_overlap)
      cross <- cross*ifelse(category[keep]==category[i],1-lambda,lambda/(length(lev)-1))
    }
    answer[i,] <- c(drop(crossprod(w,overlap%*%w)),sum(w*cross),
                    sum(abs(outer(w,w))*bound))
  }
  z <- colMeans(answer)
  c(I1=z[1],I2=z[2],score=2*z[2]-z[1],error=z[3])
}

test_that('conditional GNN handles product response integrals and finite categories', {
  old <- options(np.messages=FALSE,np.tree=FALSE)
  on.exit(options(old),add=TRUE)
  x <- c(-.83,.66,-.21,.32,-.49,.93,.09)
  y <- cbind(c(-1.4,.2,.95,-.71,1.73,-.1,.57),
             c(.42,-1.1,1.6,.03,-.63,.89,-.21),
             c(-.38,1.42,-.04,-.91,.61,-1.73,.27))
  cases <- list(gaussian2=list(p=2L,family='gaussian',degree=0L),
                epan2=list(p=2L,family='epanechnikov',degree=0L),
                gaussian3=list(p=3L,family='gaussian',degree=0L),
                mixed=list(p=1L,family='gaussian',degree=2L))
  for(name in names(cases)) {
    z <- cases[[name]]; yy <- y[,seq_len(z$p),drop=FALSE]
    category <- if(name=='mixed')factor(c('a','b','a','b','a','b','a'))else NULL
    ref <- cgnn_general_tensor_reference(x,yy,4L,2L,z$degree,z$family,category)
    expect_lt(ref['error'],2.5e-10)
    yd <- as.data.frame(yy);counts <- rep(2,z$p)
    if(!is.null(category)){yd$u <- category;counts <- c(counts,.2)}
    ctl <- if(z$degree==0L)list(regtype='lc')else
      list(regtype='lp',degree=z$degree,bernstein.basis=TRUE)
    b <- do.call(npcdensbw,c(list(xdat=data.frame(x=x),ydat=yd,
      bws=c(counts,4),bwtype='generalized_nn',bwmethod='cv.ls',
      cxkertype='gaussian',cykertype=z$family,uykertype='aitchisonaitken',
      bandwidth.compute=FALSE),ctl))
    observed <- cgnn_general_evaluate(data.frame(x=x),yd,b)
    expect_true(abs(observed-ref['score']) <= 1e-9,info=name)
  }
})


test_that('projected higher Gaussian orders retain extended and affine contracts', {
  old <- options(np.messages=FALSE,np.tree=FALSE,np.extendednn=TRUE)
  on.exit(options(old),add=TRUE)

  x <- c(-.83,.66,-.21,.32,-.49,.93,.09)
  y <- c(-1.4,.2,.95,-.71,1.73,-.1,.57)
  permutation <- c(4L,1L,7L,3L,6L,2L,5L)
  for(order in c(4L,6L,8L))for(ky in c(2L,8L))for(degree in c(0L,2L)) {
    reference <- cgnn_general_reference(x,y,4L,ky,degree,'gaussian',order)
    expect_lt(reference['error'],2.5e-10)
    ctl <- if(degree==0L)list(regtype='lc')else list(regtype='lp',degree=degree)
    xx <- data.frame(x=x[permutation]); yy <- .5+2*y[permutation]
    b <- do.call(npcdensbw,c(list(xdat=xx,ydat=yy,bws=c(ky,4),
      bwtype='generalized_nn',bwmethod='cv.ls',bandwidth.compute=FALSE,
      cykerorder=order),ctl))
    value <- cgnn_general_evaluate(xx,yy,b)
    expect_true(abs(value-reference['score']/2)<=1e-9,
      info=paste(order,ky,degree))
  }
})
