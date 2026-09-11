local({
  n <- 96L
  x <- data.frame(a=qnorm((seq_len(n)-.5)/n),
    b=qnorm(((seq_len(n)*37L)%%n+.5)/n),
    u=factor(rep(letters[1:3],length.out=n)))
  at <- data.frame(a=c(-.2,0,.2),b=c(.1,-.15,.3),
                   u=factor(c("a","b","c"),levels=levels(x$u)))
  # Dense independent oracle: actual kth distances, explicit row membership,
  # manual uniform/AA kernels, and marked sample-quantile spacing faces.
  reference <- function(data, evaluation, k, lambda=NULL) {
    n <- nrow(data); p <- length(k)
    radius <- density <- matrix(0,n,p)
    for (j in seq_len(p)) {
      z <- data[[j]]; dd<-abs(outer(z,z,"-"));diag(dd)<-Inf
      ss <- apply(dd,1,sort)
      w <- min(floor((n-1)^(2/3)),floor(min(k[j],n-1-k[j])/2))
      radius[,j]<-ss[k[j],]
      density[,j]<-2*w/(n-1)/(ss[k[j]+w,]-ss[k[j]-w,])
    }
    pilot <- function(z,t) {
      a <- sum(z<=t); w<-min(floor((n-1)^(2/3)),floor(min(a-1,n-a)/2))
      if(w<1) return(NULL)
      inds <- (a-w):(a+w)
      list(ids=inds,width=z[max(inds)]-z[min(inds)])
    }
    out <- point <- numeric(nrow(evaluation))
    for(q in seq_len(nrow(evaluation))) {
      marks <- matrix(1,n,p+as.integer(!is.null(lambda)))
      for(j in seq_len(p)) marks[,j]<-(abs(evaluation[[j]][q]-data[[j]])<radius[,j])/(2*radius[,j])
      if(!is.null(lambda)) marks[,p+1]<-ifelse(data$u==evaluation$u[q],1-lambda,lambda/2)
      A <- apply(marks,1,prod); phi<-A
      for(j in seq_len(p)) {
        I <- abs(outer(data[[j]],data[[j]],"-"))<=radius[,j]
        phi <- phi + as.vector(crossprod(I,A/(radius[,j]*density[,j])))/n
        ord<-order(data[[j]]);z<-data[[j]][ord];xx<-evaluation[[j]][q]
        probs<-mean(z<=xx)+c(-1,1)*k[j]/(n-1)
        for(h in 1:2) {
          if(probs[h]<0 || probs[h]>1) next
          if(probs[h]==0 || probs[h]==1) {phi[]<-NA_real_;break}
          s<-as.numeric(quantile(z,probs[h],type=7))
          fs<-pilot(z,s);gt<-pilot(z,(s+xx)/2)
          if(is.null(fs)||is.null(gt)) {phi[]<-NA_real_;break}
          other<-if(ncol(marks)==1L) rep(1,n) else apply(marks[,-j,drop=FALSE],1,prod)
          w<-other[ord[gt$ids]];w[c(1,length(w))]<-w[c(1,length(w))]/2
          g<-sum(w)/(n*gt$width)
          f<-(length(fs$ids)-1)/(n*fs$width)
          r<-abs(xx-s)/2
          phi<-phi+(if(h==1)-1 else 1)*g/(4*r*f)*((data[[j]]<=xx)-(data[[j]]<=s))
        }
      }
      point[q]<-mean(A);out[q]<-sqrt(var(phi)/n)
    }
    list(point=point,se=out)
  }
  run <- function() {
    old<-options(np.messages=FALSE);on.exit(options(old),add=TRUE)
    for(cols in list("a",c("a","b"),c("a","b","u"))) {
      test_that(paste("uniform ANN full boundary",paste(cols,collapse="/")), {
        d<-x[cols];e<-at[cols];p<-sum(vapply(d,is.numeric,logical(1)))
        k<-c(13L,19L)[seq_len(p)];has.cat<-"u"%in%cols
        bw<-npudensbw(dat=d,bws=c(k,if(has.cat).3),bwtype="adaptive_nn",
                       ckertype="uniform",bandwidth.compute=FALSE)
        ref<-reference(d,e,k,if(has.cat).3)
        for(tree in list(FALSE,TRUE,"auto")) {
          options(np.tree=tree)
          # The pre-existing kbandwidth constructor may mention ignored order.
          fit<-suppressWarnings(npudens(bws=bw,tdat=d,edat=e,se=TRUE))
          off<-npudens(bws=bw,tdat=d,edat=e,se=FALSE)
          expect_identical(fitted(fit),fitted(off))
          expect_equal(fitted(fit),ref$point,tolerance=2e-12)
          expect_equal(se(fit),ref$se,tolerance=2e-11)
        }
        perm<-order(sin(seq_len(n)*.713))
        changed<-suppressWarnings(npudens(bws=bw,tdat=d[perm,,drop=FALSE],edat=e,se=TRUE))
        expect_equal(se(changed),ref$se,tolerance=2e-11)
        tail.rows<-e[1:2,,drop=FALSE];tail.rows[[1]]<-c(-1.3,1.3)
        tail.ref<-reference(d,tail.rows,k,if(has.cat).3)
        tail.fit<-suppressWarnings(npudens(bws=bw,tdat=d,edat=tail.rows,se=TRUE))
        expect_true(all(is.finite(se(tail.fit))))
        expect_equal(se(tail.fit),tail.ref$se,tolerance=2e-11)
        ee<-e;ee[[1]][]<-20
        msgs<-character()
        out<-withCallingHandlers(npudens(bws=bw,tdat=d,edat=ee,se=TRUE),warning=function(w){
          msgs<<-c(msgs,conditionMessage(w));invokeRestart("muffleWarning")})
        expect_true(all(is.na(se(out))))
        expect_equal(sum(grepl("standard errors are unavailable",msgs)),1L)
        expect_true(all(fitted(out)==0))
      })
    }
  }
  if(exists(".npRmpi_with_local_regression",mode="function")) .npRmpi_with_local_regression(run()) else run()
})
