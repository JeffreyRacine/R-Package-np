test_that("convolution coordinates populate every replacement product", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force=TRUE),add=TRUE)
  old <- options(np.messages=FALSE, np.largeh=FALSE, np.largelambda=FALSE)
  on.exit(options(old), add=TRUE)
  x <- data.frame(a=c(-1.1,-.7,-.2,.1,.35,.6,.95,1.3,1.8),
                  b=c(.1,.7,-.8,1.3,-.5,.3,1.7,-1.1,.9))
  w <- cbind(1,seq_len(nrow(x))/10)
  y <- cbind(sin(x$a),cos(x$b))
  for (type in c("fixed","generalized_nn","adaptive_nn"))
    for (kernel in c("gaussian","epanechnikov","uniform"))
      for (tree in c(FALSE,TRUE)) for (loo in c(FALSE,TRUE)) {
        options(np.tree=tree)
        a <- list(txdat=x, bws=if(type=="fixed") c(.37,.63) else c(4,5),
                  bwtype=type, ckertype=kernel, bandwidth.divide=TRUE,
                  leave.one.out=loo, weights=w, tydat=y)
        for (op in c("normal","derivative","integral")) {
          all <- do.call(npksum,c(a,list(operator="convolution",
                                         permutation.operator=op)))
          expect_true(all(is.finite(all$p.ksum)))
          for (j in 1:2) {
            operators <- rep("convolution",2); operators[j] <- op
            one <- do.call(npksum,c(a,list(operator=operators)))
            expect_equal(all$p.ksum[,,,j],one$ksum,tolerance=5e-12)
          }
        }
      }
})

test_that("Gaussian replacement products agree with analytic overlap calculus", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force=TRUE),add=TRUE)
  old <- options(np.messages=FALSE,np.tree=FALSE,np.largeh=FALSE)
  on.exit(options(old),add=TRUE)
  x <- data.frame(a=c(-.8,-.3,.1,.6,1.2),b=c(.4,-.2,1.1,-.7,.8))
  h <- c(.39,.57)
  delta <- lapply(x,function(z) outer(z,z,function(t,e)e-t))
  normal <- Map(function(d,h) exp(-.5*(d/h)^2)/sqrt(2*pi)/h,delta,h)
  deriv <- Map(function(k,d,h) -d*k/h^2,normal,delta,h)
  conv <- Map(function(d,h) 0.3989422803*exp(-.25*(d/h)^2)/sqrt(2)/h,delta,h)
  fit <- npksum(txdat=x,bws=h,operator="convolution",
                 permutation.operator="derivative",bandwidth.divide=TRUE)
  expect_equal(fit$p.ksum[,1],colSums(deriv[[1]]*conv[[2]]),tolerance=3e-10)
  expect_equal(fit$p.ksum[,2],colSums(conv[[1]]*deriv[[2]]),tolerance=3e-10)
})

test_that("categorical score planes inherit continuous convolution products", {
  skip_if_not(spawn_mpi_slaves(1L), "MPI pool unavailable")
  on.exit(close_mpi_slaves(force=TRUE),add=TRUE)
  old <- options(np.messages=FALSE,np.tree=FALSE)
  on.exit(options(old),add=TRUE)
  x <- data.frame(x=c(-.8,-.3,.1,.6,1.2),
                  u=factor(c("a","b","a","c","b")))
  for (ocg in c(FALSE,TRUE)) {
    flags <- if(ocg) list(compute.ocg=TRUE) else list(compute.score=TRUE)
    got <- do.call(npksum,c(list(txdat=x,bws=c(.4,.2),
                 operator=c("convolution","normal"),
                 return.kernel.weights=TRUE,
                 return.derivative.kernel.weights=TRUE),flags))
    cont <- npksum(txdat=x[1],bws=.4,operator="convolution",
                     return.kernel.weights=TRUE)$kw
    cat <- do.call(npksum,c(list(txdat=x[2],bws=.2,
                    return.kernel.weights=TRUE,
                    return.derivative.kernel.weights=TRUE),flags))$p.kw
    expect_equal(as.double(got$p.ksum),colSums(cont*cat),tolerance=2e-12)
    expect_equal(got$p.kw,cont*cat,tolerance=2e-12)
  }
})
