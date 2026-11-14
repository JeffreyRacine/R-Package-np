local({
 run<-function(){
  old<-options(np.messages=FALSE,np.largeh=FALSE);on.exit(options(old),add=TRUE)
  test_that('fixed uniform normal factors use open physical support in every row owner',{
   for(h in c(5,49,98,107))for(shift in c(0,1024)){
    delta<-h*.Machine$double.eps
    x<-data.frame(x=shift+c(-h-delta,-h,-h+delta,0,h-delta,h,h+delta))
    at<-data.frame(x=shift)
    expected<-.5*(abs(x$x-shift)<h)
    for(tree in list(FALSE,TRUE,'auto')){
     options(np.tree=tree)
     z<-npksum(txdat=x,exdat=at,bws=h,ckertype='uniform',return.kernel.weights=TRUE)
     expect_identical(as.numeric(z$kw),as.numeric(expected))
     expect_equal(as.numeric(z$ksum),sum(expected),tolerance=1e-14)
    }
   }
  })
  test_that('fixed uniform permutation factors retain physical membership',{
   h<-c(49,98);x<-data.frame(a=c(-49,0,49,1,-1),b=c(0,98,-98,1,-1));at<-data.frame(a=0,b=0)
   norm<-lapply(1:2,function(j).5*(abs(x[[j]])<h[j]))
   cdf<-lapply(1:2,function(j)pmin(pmax((1-x[[j]]/h[j])/2,0),1))
   for(tree in list(FALSE,TRUE,'auto'))for(op in list(c('normal','normal'),c('normal','integral'),c('integral','normal'))){
    options(np.tree=tree)
    ww<-lapply(1:2,function(j)if(op[j]=='normal')norm[[j]]else cdf[[j]])
    z<-npksum(txdat=x,exdat=at,bws=h,ckertype='uniform',operator=op,
      permutation.operator='normal',return.kernel.weights=TRUE,return.derivative.kernel.weights=TRUE)
    expect_equal(as.numeric(z$kw),ww[[1]]*ww[[2]],tolerance=1e-14)
    expect_equal(as.numeric(z$p.kw[,,1]),norm[[1]]*ww[[2]],tolerance=1e-14)
    expect_equal(as.numeric(z$p.kw[,,2]),ww[[1]]*norm[[2]],tolerance=1e-14)
   }
  })
 }
 if(exists('.npRmpi_with_local_regression',mode='function')) .npRmpi_with_local_regression(run()) else run()
})
