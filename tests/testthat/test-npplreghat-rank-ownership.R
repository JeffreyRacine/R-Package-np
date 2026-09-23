test_that("partial-linear hats reject unidentified residualized designs", {
  old <- options(np.messages=FALSE, npRmpi.hat.operator.fanout=FALSE)
  on.exit(options(old),add=TRUE)
  set.seed(825); n<-29L
  z<-data.frame(z=sort(runif(n))); y<-sin(5*z$z)+rnorm(n,sd=.1)
  for(regtype in c("ll","lp")) for(shape in c("reproduced","duplicate")) {
    x<-if(shape=="reproduced") data.frame(x=z$z) else
      data.frame(x=sin(6*z$z),duplicate=sin(6*z$z))
    args<-list(xdat=x,zdat=z,ydat=y,bws=matrix(.27,ncol(x)+1L,1L),
      regtype=regtype,bandwidth.compute=FALSE)
    if(regtype=="lp")args$degree<-2L
    b<-do.call(npplregbw,args)
    expect_error(npplreg(b,txdat=x,tzdat=z,tydat=y),"rank deficient after smoothing")
    for(rows in list(seq_len(n),rev(seq_len(n)))) for(output in c("apply","matrix","constraint"))
      expect_error(npplreghat(b,txdat=x[rows,,drop=FALSE],tzdat=z[rows,,drop=FALSE],
        exdat=x[1:4,,drop=FALSE],ezdat=z[1:4,,drop=FALSE],y=y[rows],output=output),
        "npplreghat: residualized linear regressors are rank deficient")
  }
})

test_that("identified partial-linear hat outputs preserve the shared solve", {
  old<-options(np.messages=FALSE,npRmpi.hat.operator.fanout=FALSE)
  on.exit(options(old),add=TRUE)
  ns<-asNamespace(getNamespaceName(environment(npplreg)))
  complete<-get(".npreghat_complete",ns)
  set.seed(826); n<-25L
  z<-data.frame(z=sort(runif(n))); x<-data.frame(x=rnorm(n),f=factor(rep(c("10","30"),length.out=n)))
  y<-sin(5*z$z)+x$x+rnorm(n,sd=.1); yy<-cbind(y,cos(seq_len(n)))
  for(regtype in c("lc","ll","lp")) {
    args<-list(xdat=x,zdat=z,ydat=y,bws=matrix(.3,3,1),regtype=regtype,bandwidth.compute=FALSE)
    if(regtype=="lp")args$degree<-2L
    b<-do.call(npplregbw,args); ex<-x[2:6,,drop=FALSE];ez<-z[2:6,,drop=FALSE]
    X<-cbind(x$x,b$bw[[3L]]$ydati$all.dlev[[1L]][as.integer(x$f)])
    RX<-RE<-NULL
    for(j in 1:2){
      RX<-cbind(RX,X[,j]-as.vector(complete(bws=b$bw[[j+1L]],txdat=z,y=X[,j],output="apply")))
      RE<-cbind(RE,X[2:6,j]-as.vector(complete(bws=b$bw[[j+1L]],txdat=z,exdat=ez,y=X[,j],output="apply")))
    }
    Ht<-complete(bws=b$bw$yzbw,txdat=z,output="matrix")
    He<-complete(bws=b$bw$yzbw,txdat=z,exdat=ez,output="matrix")
    oracle<-He+RE%*%solve(crossprod(RX),crossprod(RX,diag(n)-Ht))
    H<-npplreghat(b,txdat=x,tzdat=z,exdat=ex,ezdat=ez,output="matrix")
    expect_equal(H,oracle,tolerance=3e-12,ignore_attr=TRUE)
    expect_equal(npplreghat(b,txdat=x,tzdat=z,exdat=ex,ezdat=ez,y=yy),H%*%yy,tolerance=3e-12)
    expect_equal(npplreghat(b,txdat=x,tzdat=z,exdat=ex,ezdat=ez,y=y,output="constraint"),
      t(H)*y,tolerance=0,ignore_attr=TRUE)
    fit<-npplreg(b,txdat=x,tzdat=z,tydat=y,exdat=ex,ezdat=ez)
    expect_equal(fitted(fit),as.vector(H%*%y),tolerance=3e-12)
    training<-npplreg(b,txdat=x,tzdat=z,tydat=y,residuals=TRUE)
    pred<-predict(training,exdat=ex,ezdat=ez,se.fit=TRUE)
    expect_equal(as.numeric(pred$fit),as.vector(H%*%y),tolerance=3e-12)
    expect_equal(as.numeric(pred$se.fit),sqrt(drop(H^2%*%as.numeric(residuals(training))^2)),tolerance=3e-12)
  }
})
