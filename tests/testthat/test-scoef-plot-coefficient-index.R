test_that("smooth coefficient plots honor the selected coefficient", {
  if (exists('spawn_mpi_slaves', mode='function')) {
    spawn_mpi_slaves()
    on.exit(close_mpi_slaves(force=TRUE), add=TRUE)
  }
  withr::local_options(np.messages=FALSE,np.plot.progress=FALSE)
  ns<-asNamespace('npRmpi');match.index<-get('.np_plot_match_coef_index',ns)
  expect_identical(match.index(3.9,3L),3L) # existing conversion is unchanged
  expect_error(match.index(4L,3L),'exceeds the 3 available coefficients')
  n<-36L;a<-sin(seq_len(n));b<-cos(seq_len(n)*.7);z<-seq(.02,.98,length.out=n)
  for(q in 1:2) {
    xd<-data.frame(a=a,b=b)[seq_len(q)];zd<-data.frame(z=z);X<-cbind(1,as.matrix(xd));beta<-c(2,3,-4)[seq_len(q+1L)];y<-drop(X%*%beta)+.01*sin(2*seq_len(n))
    bw<-npscoefbw(xdat=xd,zdat=zd,ydat=y,bws=.2,regtype='lc',bandwidth.compute=FALSE)
    fit<-npscoef(bws=bw,txdat=xd,tzdat=zd,tydat=y,betas=TRUE)
    for(k in seq_len(q+1L)) {
      p<-plot(fit,neval=5L,perspective=FALSE,coef=TRUE,coef_index=k,output='data')
      expect_identical(unique(p$targets$term),c('Intercept',names(xd))[k])
      for(i in seq_along(p$panels)) {
        vars<-p$panels[[i]]$variables;zc<-vars$coordinate[vars$name=='z'];e<-as.data.frame(p,panel=i)
        expect_length(zc,1L);expect_equal(nrow(e),5L)
        oracle<-vapply(e[[zc]],function(at){w<-exp(-.5*((z-at)/.2)^2);solve(crossprod(X,X*w),crossprod(X,y*w))[k]},0)
        expect_length(oracle,nrow(e));expect_true(all(is.finite(oracle)))
        expect_equal(e$estimate,oracle,tolerance=1e-10)
      }
    }
    expect_error(plot(fit,neval=5L,coef=TRUE,coef_index=q+2L,output='data'),'available coefficients')
    expect_error(plot(bw,xdat=xd,ydat=y,zdat=zd,neval=5L,coef=TRUE,coef_index=q+2L,output='data'),'available coefficients')
    mean1<-plot(fit,neval=5L,perspective=FALSE,coef_index=1L,output='data')
    mean99<-plot(fit,neval=5L,perspective=FALSE,coef_index=99L,output='data')
    for(i in seq_along(mean1$panels))expect_identical(as.data.frame(mean1,panel=i),as.data.frame(mean99,panel=i))
  }
})
