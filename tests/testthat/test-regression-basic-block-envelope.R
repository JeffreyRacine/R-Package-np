# A literal count-law oracle independent of inference reducers/studentizers.
test_that("basic block envelopes preserve selected scale and literal counts", {
  old <- options(np.messages=FALSE);on.exit(options(old),add=TRUE)
  set.seed(98121)
  configs <- list(
    ordinary=data.frame(x=runif(120)),
    q2=data.frame(x=runif(160),z=runif(160)),
    binary=data.frame(x=runif(160),u=factor(rep(1:2,80))),
    pooled=data.frame(x=runif(180),u=factor(rep(1:3,60))),
    categorical=data.frame(u=factor(rep(1:3,30))))
  variance <- function(F,L,law) {
    F<-F-rowMeans(F);n<-ncol(F)
    if(L==1)return(rowSums(F^2))
    if(law=="fixed") {
      sums<-function(w)if(w==0)rep(0,nrow(F)) else {
        Z<-sapply(seq_len(n),function(i)rowSums(F[,((i-1L+seq_len(w)-1L)%%n)+1L,drop=FALSE]));rowMeans(Z^2)
      }
      return((n%/%L)*sums(L)+sums(n%%L))
    }
    v<-rowSums(F^2)
    for(k in seq_len(n-1L))v<-v+2*(n-k)/n*(1-1/L)^k*rowSums(F*F[,((seq_len(n)-1L+k)%%n)+1L,drop=FALSE])
    v
  }
  bands <- function(R,c,alpha) {
    B<-ncol(R);m<-nrow(R);rs<-t(apply(R,1,sort));aa<-t(apply(abs(R),1,sort))
    out<-setNames(lapply(1:3,function(i)array(0,c(m,2,length(alpha)))),c('pointwise','simultaneous','bonferroni'))
    for(a in seq_along(alpha)) {
      k<-1+floor(B*alpha[a]/2);kb<-1+floor(B*alpha[a]/(2*m));ks<-ceiling(B*(1-alpha[a]))
      S<-c(sort(apply(R,2,min))[k],sort(apply(R,2,max))[B+1-k]);q<-sort(apply(abs(R),2,max))[ks]
      for(b in names(out)) {
        lim<-switch(b,pointwise=rs[,c(k,B+1-k)],simultaneous=cbind(rep(S[1],m),rep(S[2],m)),bonferroni=cbind(pmin(rs[,kb],S[1]),pmax(rs[,B+1-kb],S[2])))
        ref<-switch(b,pointwise=aa[,ks],simultaneous=rep(q,m),bonferroni=pmax(aa[,ceiling(B*(1-alpha[a]/m))],q))
        out[[b]][,,a]<-cbind(c-pmax(lim[,2],ref),c-pmin(lim[,1],-ref))
      }
    };out
  }
  for(kind in names(configs)) {
    x<-configs[[kind]];n<-nrow(x);nc<-which(vapply(x,is.numeric,logical(1)))
    y<-if(length(nc))sin(3*x[[nc[1]]]) else as.numeric(x[[1]])
    y<-y+rnorm(n,sd=.4)
    bw<-npregbw(xdat=x,ydat=y,bws=rep(.3,ncol(x)),regtype=if(kind=="categorical")"lc"else"ll",bandwidth.compute=FALSE)
    grid<-x[seq_len(if(kind=='ordinary')26L else 7L),,drop=FALSE]
    targets <- c('mean', if(length(nc))paste0('d',seq_along(nc)),
      if(length(nc)<ncol(x))paste0('c',which(!vapply(x,is.numeric,logical(1)))))
    for(target in targets)for(law in c('fixed','geom'))for(L in c(1L,5L)) {
      state<-.np_bias_regression_prepare(bw,x,y,np_bias_control(),adapt=FALSE,law=law)
      C<-do.call(rbind,lapply(split(seq_len(nrow(grid)),ceiling(seq_len(nrow(grid))/25)),function(ii)
        .np_bias_regression_rows(state,grid[ii,,drop=FALSE],data.frame(row=seq_along(ii),target=target))$W))
      v<-.np_bias_regression_rows(state,grid[1:min(25,nrow(grid)),,drop=FALSE],data.frame(row=1:min(25,nrow(grid)),target=target))
      se0<-v$se
      if(nrow(grid)>25)se0<-c(se0,.np_bias_regression_rows(state,grid[26,,drop=FALSE],data.frame(row=1L,target=target))$se)
      F<-sweep(C,2,state$residual,'*');F<-F-rowMeans(F)
      fit<-if(!is.null(state$main.capture) || !is.null(state$local.inference))
        .np_plot_regression_eval(bw,x,y,grid,gradients=target!="mean",need.bias.capture=TRUE) else NULL
      B<-199L;alpha<-c(.05,.2)
      z<-.np_with_seed(98122,.np_bias_inference(state,grid,targets=target,B=B,law=law,fit=fit,blocklen=L,calibration='basic-envelope',refinement='pilot-block',alpha=alpha,resolution.warn=FALSE))
      d<-.np_bias_stream(z$rng$start,.np_block_counts_drawer(n,B,L,law)(1L,B)-1)
      R<-F%*%d$value
      if(isTRUE(state$pooled)) {
        expect_equal(z$se,sqrt(variance(F,L,law)),tolerance=1e-10)
        expect_null(z$studentization)
        expect_match(z$covariance$studentization,'observed')
      } else {
        v1 <- variance(F,1L,law)
        expect_true(all(se0[v1==0]==0))
        anchored<-se0*sqrt(variance(F,L,law)/ifelse(v1==0,1,v1))
        expect_equal(z$se,anchored,tolerance=1e-10)
        denominator <- apply(R,1,sd)
        expect_true(all(se0[denominator==0]==0))
        R<-anchored*(R/ifelse(denominator==0,1,denominator))
      }
      if(isTRUE(state$pooled)) {
        expect_identical(z$root.units,"endpoint units; raw centered pooled block scores")
        expect_identical(z$provider.refinement$owner,"exact conditional pooled block score variance")
        if(L==1L) {
          retained<-.np_with_seed(98122,.np_bias_inference(state,grid,targets=target,B=B,law=law,fit=fit,blocklen=L,
            calibration='tail-envelope',refinement='pilot-block',resolution.warn=FALSE))
          expect_identical(retained$studentization$construction,"conditional-block-score-v1")
        }
      } else {
        expect_identical(z$root.units,"endpoint units; pilot-anchored block residual scores")
        expect_identical(z$provider.refinement$owner,"pilot times selected-block to iid covariance ratio")
      }
      expect_identical(z$rng$end,d$end)
      expected<-bands(R,z$center,alpha)
      for(b in names(expected))for(k in seq_along(alpha))
        expected[[b]][z$se==0,,k] <- z$center[z$se==0]
      for(b in names(expected))expect_equal(unname(z$bands[[b]]),unname(expected[[b]]),tolerance=1e-10,info=paste(kind,target,law,L,b))
      expect_true(all(z$bands$simultaneous[,1,]<=z$bands$pointwise[,1,]))
      expect_true(all(z$bands$simultaneous[,2,]>=z$bands$pointwise[,2,]))
    }
    fixed.state<-.np_bias_regression_prepare(bw,x,y,np_bias_control(),adapt=FALSE,law='fixed')
    expect_error(.np_with_seed(98122,.np_bias_inference(fixed.state,grid,B=99L,law='fixed',fit=fit,blocklen=5L,calibration='basic-envelope',resolution.warn=FALSE)), "requires refinement='pilot-block'")
    expect_error(.np_with_seed(98122,.np_bias_inference(state,grid,B=1L,law='geom',fit=fit,blocklen=5L,calibration='basic-envelope',refinement='pilot-block',resolution.warn=FALSE)), "at least two draws")
    expect_error(.np_with_seed(98122,.np_bias_inference(fixed.state,grid,B=99L,law='fixed',fit=fit,blocklen=n,calibration='basic-envelope',refinement='pilot-block',resolution.warn=FALSE)), 'nonpositive|undefined')
  }
})


test_that("public basic block admission preserves state and panel endpoints", {
  old<-options(np.messages=FALSE,np.plot.progress=FALSE);on.exit(options(old),add=TRUE)
  set.seed(98125);x<-data.frame(x=runif(90));y<-sin(3*x$x)+rnorm(90,sd=.3)
  bw<-npregbw(xdat=x,ydat=y,bws=.3,regtype="ll",bandwidth.compute=FALSE)
  fit<-npreg(bws=bw,txdat=x,tydat=y);held<-fit;rng<-.Random.seed
  for(law in c("fixed","geom")) {
    out<-suppressWarnings(plot(fit,center="bias-corrected",errors="bootstrap",bootstrap=law,
      B=199L,neval=4L,band="all",output="data",random.seed=98126L,
      boot.control=np_boot_control(calibration="basic-envelope",blocklen=3L)))
    z<-as.data.frame(out,band="simultaneous")
    expect_equal(nrow(z),4L);expect_true(all(is.finite(z$lower)&is.finite(z$upper)))
    expect_true(all(z$lower<=z$upper));expect_identical(fit,held);expect_identical(.Random.seed,rng)
  }
})

test_that("geometric score variance matches enumerated index paths", {
  F<-rbind(c(-1,0,1),c(1,-2,1));n<-ncol(F);paths<-as.matrix(expand.grid(rep(list(seq_len(n)),n)))
  for(L in c(1L,2L,3L)) {
    weights<-rep(1/n,nrow(paths))
    for(k in 2:n)weights<-weights*((1-1/L)*(paths[,k]==paths[,k-1L]%%n+1L)+1/(L*n))
    sums<-sapply(seq_len(nrow(paths)),function(i)rowSums(F[,paths[i,],drop=FALSE]))
    exact<-drop((sums^2)%*%weights)
    expect_equal(sum(weights),1,tolerance=1e-14)
    expect_equal(.np_bias_score_variance(F,L,"geom"),exact,tolerance=1e-12)
  }
})
