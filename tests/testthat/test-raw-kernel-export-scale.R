test_that("raw kernel exports are independent of sum normalization and power", {
  old <- options(np.messages=FALSE,np.largeh=FALSE,np.largelambda=FALSE)
  on.exit(options(old),add=TRUE)
  x <- data.frame(a=c(-1.1,-.7,-.2,.1,.35,.6,.95,1.3,1.8),
                  b=c(.1,.7,-.8,1.3,-.5,.3,1.7,-1.1,.9))
  for(type in c("fixed","generalized_nn","adaptive_nn"))
    for(tree in c(FALSE,TRUE)) for(loo in c(FALSE,TRUE)) {
      options(np.tree=tree)
      a <- list(txdat=x,bws=if(type=="fixed")c(.37,.63)else c(4,5),
                bwtype=type,leave.one.out=loo,
                return.kernel.weights=TRUE,return.derivative.kernel.weights=TRUE)
      for(op in c("normal","derivative","integral","convolution")) {
        raw <- do.call(npksum,c(a,list(operator=op,
          permutation.operator="derivative",bandwidth.divide=FALSE)))
        for(power in c(1,2)) {
          got <- do.call(npksum,c(a,list(operator=op,
            permutation.operator="derivative",bandwidth.divide=TRUE,kernel.pow=power)))
          expect_equal(got$kw,raw$kw,tolerance=2e-12)
          expect_equal(got$p.kw,raw$p.kw,tolerance=2e-12)
          noexport <- a
          noexport$return.kernel.weights <- FALSE
          noexport$return.derivative.kernel.weights <- FALSE
          sums <- do.call(npksum,c(noexport,list(operator=op,
            permutation.operator="derivative",bandwidth.divide=TRUE,kernel.pow=power)))
          expect_identical(got$ksum,sums$ksum)
          expect_identical(got$p.ksum,sums$p.ksum)
        }
      }
    }
})

test_that("raw convolution exports have the physical Gaussian overlap scale", {
  old <- options(np.messages=FALSE,np.tree=FALSE,np.largeh=FALSE)
  on.exit(options(old),add=TRUE)
  x <- c(-.8,-.3,.1,.6,1.2)
  for(type in c("fixed","generalized_nn","adaptive_nn")) {
    he <- if(type=="fixed")rep(.39,length(x))else vapply(x,function(z)sort(abs(x-z))[4],0)
    ht <- if(type=="adaptive_nn")vapply(x,function(z)sort(abs(x-z))[5],0)else he
    expected <- outer(seq_along(x),seq_along(x),Vectorize(function(i,j)
      0.3989422803*ht[i]*he[j]/sqrt(ht[i]^2+he[j]^2)*
        exp(-.5*(x[i]-x[j])^2/(ht[i]^2+he[j]^2))))
    for(divide in c(FALSE,TRUE)) {
      got <- npksum(txdat=data.frame(x),exdat=data.frame(x),
        bws=if(type=="fixed").39 else 4,bwtype=type,operator="convolution",
        bandwidth.divide=divide,return.kernel.weights=TRUE)
      expect_equal(got$kw,expected,tolerance=3e-12)
    }
  }
})

test_that("categorical replacement exports retain raw convolution scale", {
  old <- options(np.messages=FALSE,np.tree=FALSE,np.largeh=FALSE)
  on.exit(options(old),add=TRUE)
  x <- data.frame(x=c(-.8,-.3,.1,.6,1.2),u=factor(c("a","b","a","c","b")))
  for(type in c("fixed","generalized_nn","adaptive_nn"))
    for(ocg in c(FALSE,TRUE)) {
      a <- c(list(txdat=x,bws=if(type=="fixed")c(.4,.2)else c(3,.2),
        bwtype=type,operator=c("convolution","normal"),
        return.kernel.weights=TRUE,return.derivative.kernel.weights=TRUE),
        if(ocg)list(compute.ocg=TRUE)else list(compute.score=TRUE))
      raw <- do.call(npksum,c(a,list(bandwidth.divide=FALSE)))
      normalized <- do.call(npksum,c(a,list(bandwidth.divide=TRUE)))
      expect_equal(normalized$kw,raw$kw,tolerance=2e-12)
      expect_equal(normalized$p.kw,raw$p.kw,tolerance=2e-12)
    }
})
