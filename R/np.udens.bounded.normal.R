# Complete normal-scale inference for fixed-bound q1 density plots only.
# The original guide owns the center; the ratio-pilot guide owns the variance.
.np_udens_bounded_normal_request <- function(bws, method, band, calibration, alpha) {
  if (!identical(band, "pmzsd")) return(FALSE)
  if (!identical(method, "asymptotic") || !isTRUE(bws[["ndim", exact=TRUE]] == 1L))
    stop("bounded density pmzsd requires q=1 and errors='asymptotic'", call.=FALSE)
  if (!is.null(calibration))
    stop("bounded density asymptotic pmzsd does not accept bootstrap calibration", call.=FALSE)
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= .5)
    stop("normal bias intervals require alpha in (0, 0.5)", call.=FALSE)
  TRUE
}

.np_udens_bounded_normal_analytic <- function(theta) {
  a <- theta[1L]+1; b <- theta[2L]+1
  if (any(!is.finite(c(a,b))) || a <= 0 || b <= 0)
    stop("invalid bounded normal reference guide", call.=FALSE)
  s <- a+b; u <- trigamma(s)
  list(A=lbeta(a,b), mean=c(digamma(a)-digamma(s),digamma(b)-digamma(s)),
    I=matrix(c(trigamma(a)-u,-u,-u,trigamma(b)-u),2L,2L))
}

.np_udens_bounded_normal_fit <- function(target, theta) {
  for (it in 0:79) {
    st <- .np_udens_bounded_normal_analytic(theta)
    grad <- st$mean-target
    ev <- eigen(st$I,symmetric=TRUE,only.values=TRUE)$values
    if (any(!is.finite(ev)) || min(ev) <= 0 || max(ev)/min(ev) > 1e13)
      stop("bounded normal reference information conditioning", call.=FALSE)
    d <- solve(st$I,grad); dec <- sum(grad*d); obj <- st$A-sum(theta*target)
    # Also satisfy the caller's existing relative score certificate. An
    # absolute residual below 2e-8 alone can stop too soon for small targets.
    if (dec < 1e-14 && max(abs(grad)) < 2e-8 &&
        max(abs(grad)) <= 1e-12 + 2e-8 * max(abs(target))) break
    step <- 1
    for (j in 1:2) if (d[j] > 0) step <- min(step,.99*(theta[j]+1)/d[j])
    accepted <- FALSE
    for (ls in 1:40) {
      trial <- theta-step*d; s <- .np_udens_bounded_normal_analytic(trial)
      if (s$A-sum(trial*target) <= obj-1e-4*step*dec+1e-14) {
        theta <- trial; accepted <- TRUE; break
      }
      step <- step*.5
    }
    if (!accepted) stop("bounded normal reference score line search", call.=FALSE)
    if (it == 79L) stop("bounded normal reference score stationarity", call.=FALSE)
  }
  list(theta=theta,iterations=it,
    score_error=max(abs(.np_udens_bounded_normal_analytic(theta)$mean-target)))
}

.np_udens_bounded_normal_budget <- function(n, nodes, targets) {
  # Account for simultaneous input, weighted-score and two-level output copies.
  bytes <- 8*(12*32*max(n,nodes)+12*(n+nodes)+30*targets)
  if (!is.finite(bytes) || bytes > 192*1024^2)
    stop("bounded normal plot workspace exceeds 192 MiB", call.=FALSE)
  bytes
}

.np_udens_bounded_normal_pilot <- function(x, theta, outer, targets) {
  st <- .np_udens_beta_marginal_rule(theta,outer)
  budget <- .np_udens_bounded_normal_budget(length(x),st$nodes,targets)
  ex <- .np_udens_bounded_normal_analytic(theta)
  .np_udens_beta_check(st$mean,ex$mean,"analytic beta mean",1e-10)
  .np_udens_beta_check(st$I,ex$I,"analytic beta information")
  Td <- cbind(log(x),log1p(-x)); lgx <- drop(Td%*%theta)-st$A
  n <- length(x); rho <- n^(-2/5); lp <- numeric(st$nodes)
  for (start in seq.int(1L,st$nodes,by=32L)) {
    ix <- seq.int(start,min(st$nodes,start+31L)); v <- exp(st$lv[ix])
    a <- v/rho; b <- (1-v)/rho
    logs <- outer(a,log(x))+outer(b,log1p(-x))-lbeta(a+1,b+1)
    logs <- sweep(logs,2,lgx)
    lp[ix] <- apply(logs,1,.np_udens_beta_logsum)-log(n)
  }
  lw <- log(st$prob)+lp; normal <- .np_udens_beta_logsum(lw)
  pp <- exp(lw-normal); pp <- pp/sum(pp)
  if (!is.finite(exp(normal)) || exp(normal) <= 0 || any(!is.finite(pp)))
    stop("invalid bounded normal pilot normalization", call.=FALSE)
  target <- drop(crossprod(pp,st$T))
  fit <- .np_udens_bounded_normal_fit(target,theta)
  sq <- .np_udens_beta_marginal_rule(fit$theta,outer)
  budget <- max(budget,.np_udens_bounded_normal_budget(n,sq$nodes,targets))
  .np_udens_beta_check(sq$mean,target,"reference score equation",2e-8)
  list(st=st,prob=pp,normal=exp(normal),target=target,fit=fit,sq=sq,
    V=crossprod(sweep(st$T,2,target),sweep(st$T,2,target)*pp),workspace=budget)
}

.np_udens_bounded_normal_variance <- function(A, T, prob, D, I, n) {
  muA <- drop(A%*%prob); muT <- drop(crossprod(prob,T))
  Ac <- A-muA; Tc <- sweep(T,2,muT)
  extra <- (D%*%solve(I))%*%t(Tc); psi <- Ac+extra
  value <- rowSums(sweep(psi^2,2,prob,'*'))/n
  VA <- rowSums(sweep(Ac^2,2,prob,'*'))/n
  VT <- rowSums(sweep(extra^2,2,prob,'*'))/n
  cross <- 2*rowSums(sweep(Ac*extra,2,prob,'*'))/n
  check <- .np_udens_beta_check(value,VA+VT+cross,"expanded full variance")
  if (any(!is.finite(value)) || any(value < 0))
    stop("invalid bounded normal complete variance", call.=FALSE)
  list(value=value,meanA=muA,check=check)
}

.np_udens_bounded_normal_moments <- function(x,h,grid,kernel,
                                            progress=function(...) NULL) {
  if (!is.numeric(x) || length(x)<10L || !is.numeric(grid) || !length(grid) ||
      any(!is.finite(c(x,grid,h))) || any(x<=0|x>=1) || any(grid<=0|grid>=1) ||
      length(h)!=1L || h<=0 || length(kernel)!=1L || !kernel%in%c('beta','gaussian'))
    stop("invalid bounded normal q1 moment request", call.=FALSE)
  n <- length(x)
  .np_udens_bounded_normal_budget(n,0,length(grid))
  progress('original guide',0L,1L)
  theta <- .np_udens_beta_q1_fit(x)$theta
  progress('original guide',1L,1L)
  Tg <- cbind(log(grid),log1p(-grid))
  score <- if (kernel=='beta') .np_udens_beta_q1_score else .np_udens_fgd_score
  inner <- if (kernel=='beta') c(64L,128L) else c(24L,48L)
  outer <- c(128L,256L); result <- vector('list',2L)
  variance.plan <- .np_udens_q1_operator_plan(theta,h,grid)
  for (level in 1:2) {
    progress(paste('reference pilot',level),0L,1L)
    p <- .np_udens_bounded_normal_pilot(x,theta,outer[level],length(grid))
    st <- p$st; sq <- p$sq; tq <- p$fit$theta
    progress(paste('reference pilot',level),1L,1L)
    vr <- .np_udens_q1_operator_rule(theta,variance.plan,outer[level])
    qp <- .np_udens_q1_operator_plan(tq,h,grid)
    qr <- .np_udens_q1_operator_rule(tq,qp,outer[level],keep.log=TRUE)
    for(key in c("A","mean","I")) .np_udens_beta_check(qr[[key]],sq[[key]],
      paste("normal operator variance guide",key),if(key=="mean")1e-10 else 1e-9)
    extra.workspace <- .np_udens_q1_variance_workspace(n,min(32L,length(grid)),
      qr$nodes,vr$nodes)
    if(extra.workspace>192*1024^2)stop('bounded normal plot workspace exceeds 192 MiB',call.=FALSE)
    p$workspace <- max(p$workspace,extra.workspace)
    for(key in c('A','mean','I')) .np_udens_beta_check(vr[[key]],st[[key]],
      paste('normal variance guide',key),if(key=='mean')1e-10 else 1e-9)
    vw <- .np_udens_q1_ratio_weights(x,theta,st$A,vr)
    .np_udens_beta_check(vw$normal,p$normal,'normal cross-rule normalization')
    # The variance-guide nodes provide a second evaluation of the SAME Q law.
    qst<-qr
    qst$prob<-exp(qr$log.prob+drop(qr$T%*%(theta-tq))+qr$A-st$A)
    .np_udens_beta_check(sum(qst$prob),1,'normal cross-guide mass',1e-10)
    qw<-.np_udens_q1_ratio_weights(x,theta,st$A,qst)
    .np_udens_beta_check(qw$normal,p$normal,'normal cross-guide normalization')
    fitted <- exp(drop(Tg%*%theta)-st$A)
    fittedQ <- exp(drop(Tg%*%tq)-sq$A)
    Gq <- sweep(sweep(Tg,2,sq$mean),1,fittedQ,'*')
    fields <- c('center','se','mu','fitted','meanQ')
    z <- setNames(lapply(fields,function(i) numeric(length(grid))),fields)
    z$theta <- theta; z$thetaQ <- tq
    z$reference <- list(normal=p$normal,score_error=p$fit$score_error,
      iterations=p$fit$iterations,V=p$V,I=sq$I)
    z$workspace <- p$workspace; z$variance.checks <- list()
    for (first in seq.int(1L,length(grid),by=32L)) {
      ix <- seq.int(first,min(length(grid),first+31L)); g <- grid[ix]
      progress(paste('complete moments',level),first-1L,length(grid))
      coordinates <- list(vr,qr,list(lv=log(x),l1v=log1p(-x)))
      columns <- .np_udens_q1_score_columns(g,coordinates,h,n^(-7/17),inner[level],score)
      A <- columns[[1L]]; Aq <- columns[[2L]]; obs <- columns[[3L]]
      Av <- A
      rm(columns,coordinates)
      mu <- drop(A%*%vr$prob)
      Bq <- Aq%*%(sweep(qr$T,2,sq$mean)*qr$prob)
      Dq <- Gq[ix,,drop=FALSE]-Bq
      original.mean <- drop(Aq%*%qw$prob)
      vq <- .np_udens_bounded_normal_variance(Av,vr$T,vw$prob,Dq,sq$I,n)
      .np_udens_beta_check(vq$meanA,original.mean,'normal cross-rule mean')
      z$center[ix] <- rowMeans(obs)-mu+fitted[ix]
      z$se[ix] <- sqrt(vq$value); z$mu[ix] <- mu
      z$fitted[ix] <- fitted[ix]; z$meanQ[ix] <- original.mean
      z$variance.checks[[length(z$variance.checks)+1L]] <- vq$check
      progress(paste('complete moments',level),max(ix),length(grid))
    }
    result[[level]] <- z
  }
  fields <- c('center','se','mu','fitted','meanQ','thetaQ')
  checks <- lapply(fields,function(name) .np_udens_beta_check(
    result[[1L]][[name]],result[[2L]][[name]],paste('normal refinement',name),
    if (name=='mu') 1e-10 else 1e-9))
  names(checks) <- fields
  result[[1L]]$checks <- checks
  result[[1L]]
}

.np_udens_bounded_normal_public <- function(state,grid,alpha,band,kernel) {
  if (!identical(band,'pmzsd') || !identical(state$q,1L) ||
      !identical(state$variant,'asymptotic'))
    stop("invalid bounded normal plot state",call.=FALSE)
  if (!is.data.frame(grid) || !nrow(grid) || !identical(names(grid),names(state$x)) ||
      !all(vapply(grid,is.numeric,logical(1L))) || any(!is.finite(as.matrix(grid))))
    stop("invalid bounded normal evaluation frame",call.=FALSE)
  normalized <- (grid[[1L]]-state$lower)/state$width
  if (any(normalized<=0|normalized>=1))
    stop("bounded normal bands require strictly interior evaluation points; literal support endpoints are excluded",call.=FALSE)
  progress <- function(stage,done,total) .np_plot_progress_notify(
    stage=paste('Bounded density normal:',stage),done=done,total=total,unit='work')
  z <- .np_plot_corrected_activity_run('Constructing fixed-bound density normal intervals',
    .np_udens_bounded_normal_moments(state$normalized[,1L],state$h,normalized,kernel,progress))
  inference <- list(center=z$center/state$jacobian,se=z$se/state$jacobian,bands=list(),
    settings=list(construction='density-fixed-q1-ratio-pilot-normal-v1',
      calibration='normal-quantile',law='asymptotic',kernel=kernel,q=1L,n=state$n,
      alpha=alpha,B=0L,available.bands='pmzsd',grid.size=nrow(grid),
      family='marginal intervals at supplied interior targets',
      bounds=cbind(state$lower,state$lower+state$width),bandwidth=state$h*state$width,
      calibration.status='complete ratio-pilot variance; unchanged corrected center'),
    reference=list(original.guide=z$theta,variance.guide=z$thetaQ,
      normalization=z$reference$normal,score.error=z$reference$score_error,
      score.iterations=z$reference$iterations,score.covariance=z$reference$V,
      guide.information=z$reference$I),
    numerical.checks=list(refinement=z$checks,variance=z$variance.checks,workspace=z$workspace))
  inference <- .np_plot_bias_normal(inference,alpha)
  inference$intervals <- inference$absolute.intervals <- inference$bands
  endpoints <- inference$bands$pmzsd
  list(boot.err=cbind(inference$center-endpoints[,1L,1L],
    endpoints[,2L,1L]-inference$center,inference$center),
    boot.all.err=NULL,inference=inference,bxp=list())
}
