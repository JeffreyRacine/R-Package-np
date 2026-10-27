# Private F4 positive marginal covariance correction. The estimating center,
# empirical replicate pivots and correlation calibration owner are unchanged.
.np_udens_b4_pilot_moment <- function(state,left,right=NULL,tau=2,corrected=TRUE) {
 numeric<-vapply(state$x,is.numeric,logical(1));q<-sum(numeric)
 x<-state$x[numeric];h<-state$h[numeric];t<-as.double(left[numeric]);weight<-rep(1,state$n)
 if(is.null(right)){
  radius<-numeric(state$n)
  for(j in seq_len(q)){u<-(t[j]-x[[j]])/h[j];weight<-weight*dnorm(u/sqrt(1+tau))/(sqrt(1+tau)*h[j]);radius<-radius+u^2}
  value<-weight*if(corrected)(1+q/(2*(1+tau))-radius/(2*(1+tau)^2))else 1
 }else{
  s<-as.double(right[numeric]);et<-es<-cross<-numeric(state$n);v<-tau/(1+2*tau)
  for(j in seq_len(q)){
   a<-(t[j]-x[[j]])/h[j];b<-(s[j]-x[[j]])/h[j]
   weight<-weight*exp(-(a-b)^2/4-(a+b)^2/(4*(1+2*tau)))/(2*pi*h[j]^2*sqrt(1+2*tau))
   mt<-((1+tau)*a-tau*b)/(1+2*tau);ms<-((1+tau)*b-tau*a)/(1+2*tau)
   et<-et+mt^2+v;es<-es+ms^2+v;cross<-cross+4*mt*ms*v+2*v^2
  }
  value<-weight*if(corrected)((1+q/2-et/2)*(1+q/2-es/2)+cross/4)else 1
 }
 if(state$r+state$s>0)for(z in state$category){
  target<-as.double(toMatrix(left[z$index]));if(!is.null(right))target<-c(target,as.double(toMatrix(right[z$index])))
  value<-value*.np_udens_b2_category_moment(z,target,corrected)
 }
 if(any(!is.finite(value)))stop('nonfinite F4 pilot moment',call.=FALSE)
 value
}
.np_udens_b4_covariance <- function(state,grid,model,corrected=TRUE) {
 if(!identical(state[['metadata',exact=TRUE]][['type',exact=TRUE]],'fixed'))stop('F4 requires fixed bandwidth',call.=FALSE)
 .np_udens_bc_grid(state,grid)
 q<-sum(vapply(state$x,is.numeric,logical(1)))
 if(q==0L)return(list(covariance=model$covariance,se=model$se,ratio=rep(1,nrow(grid)),wide.variance=model$se^2))
 m<-nrow(grid);.np_udens_bc_budget(state,m,1L,25L,TRUE);wide<-numeric(m)
 for(i in seq_len(m)){
  a<-.np_udens_b4_pilot_moment(state,grid[i,,drop=FALSE],tau=2,corrected=corrected)
  b<-.np_udens_b4_pilot_moment(state,grid[i,,drop=FALSE],grid[i,,drop=FALSE],tau=2,corrected=corrected)
  wide[i]<-(mean(b)-mean(a)^2)/state$n
  .np_plot_progress_notify(stage = "Correcting density scales",
    done = i, total = m, unit = "targets")
 }
 if(any(!is.finite(wide))||any(wide<=0))stop('F4 requires positive finite widened-pilot variances',call.=FALSE)
 d<-sqrt(diag(model$covariance)/wide)
 if(any(!is.finite(d))||any(d<=0))stop('F4 requires positive finite scale factors',call.=FALSE)
 C<-model$covariance*tcrossprod(d);se<-sqrt(diag(C));ratio<-se/model$se
 if(any(!is.finite(C))||any(!is.finite(ratio))||any(ratio<=0))stop('invalid F4 corrected covariance',call.=FALSE)
 list(covariance=C,se=se,ratio=ratio,wide.variance=wide)
}
.np_udens_b4_adjust <- function(state,grid,base,corrected=TRUE) {
 if(sum(vapply(state$x,is.numeric,logical(1)))==0L)return(base)
 model<-.np_udens_b4_covariance(state,grid,base,corrected)
 map<-function(ints)lapply(ints,function(z){out<-sweep(sweep(sweep(z,1L,base$center,'-'),1L,model$ratio,'*'),1L,base$center,'+');if(anyNA(out)||any(out[,1,]>out[,2,]))stop('invalid F4 endpoint ordering',call.=FALSE);out})
 out<-base;out$intervals<-map(base$intervals);out$absolute.intervals<-map(base$absolute.intervals);out$covariance<-model$covariance;out$se<-model$se
 out$settings$construction<-'density-bc-v4-positive-marginal-covariance';out$settings$observed.scale<-'positive marginal Richardson congruence; continuous pilot tau=2';out$settings$correlation<-'unchanged F3 positive-pilot correlation';out$F4<-list(ratio=model$ratio,wide.variance=model$wide.variance)
 out
}

# Keep the complete operation in the existing transactional RNG contract.
.np_udens_b4_inference <- function(state, grid, law = c("gaussian", "inid"),
                                    alpha = c(.05, .1), B = 9999L,
                                    tile = 25L, corrected = TRUE) {
  existed <- exists(".Random.seed", .GlobalEnv, inherits = FALSE)
  previous <- if (existed) get(".Random.seed", .GlobalEnv, inherits = FALSE) else NULL
  kind <- RNGkind(); success <- FALSE
  on.exit({
    if (!success) {
      do.call(RNGkind, as.list(kind))
      if (existed) assign(".Random.seed", previous, .GlobalEnv)
      else if (exists(".Random.seed", .GlobalEnv, inherits = FALSE))
        rm(list = ".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  law <- match.arg(law)
  base <- if (law == "inid")
    .np_udens_b3_inid(state, grid, alpha, B, tile, corrected) else
    .np_udens_b2_gaussian(state, grid, alpha, B, tile, corrected)
  out <- .np_udens_b4_adjust(state, grid, base, corrected)
  if (law == "gaussian") out$absolute.intervals <- NULL
  out$settings$construction <- "density-bc-v4-positive-marginal-covariance"
  out$settings$observed.scale <- "positive marginal Richardson congruence; continuous pilot tau=2"
  out$settings$correlation <- "unchanged F3 positive-pilot correlation"
  success <- TRUE
  out
}
