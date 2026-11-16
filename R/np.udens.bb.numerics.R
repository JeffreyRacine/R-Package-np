# Numerical owner for the explicit plot-only beta-binomial bootstrap law.
# Existing U5/Q6/beta paths do not call these helpers.
.np_udens_bb_check <- function(a, b, label) {
  if (!identical(dim(a),dim(b)) || length(a)!=length(b) ||
      !length(a) || any(!is.finite(c(a,b))))
    stop(paste("invalid beta-binomial comparison:",label),call.=FALSE)
  ratio <- max(abs(a-b)/(1e-12+1e-9*abs(b)))
  if (!is.finite(ratio) || ratio>1)
    stop(paste("beta-binomial refinement failed:",label,ratio),call.=FALSE)
  ratio
}

.np_udens_bb_operator <- function(grid,h) {
  if (!is.numeric(grid) || !length(grid) || any(!is.finite(grid)) ||
      any(grid<=0|grid>=1) || !is.numeric(h) || length(h)!=1L || !is.finite(h) || h<=0)
    stop("beta-binomial operator requires interior targets and a positive bandwidth",call.=FALSE)
  rule <- .np_udens_q6_rule(16L)
  t(vapply(grid,function(g) {
    l <- max(-sqrt(5),-g/h);r <- min(sqrt(5),(1-g)/h)
    half <- (r-l)/2;mid <- (r+l)/2;u <- mid+half*rule$x
    P <- .np_udens_u5_poly(rule$x,5L)
    w <- half*rule$w*(3/(4*sqrt(5)))*(1-u^2/5)
    G <- crossprod(P,P*w);p0 <- drop(.np_udens_u5_poly(-mid/half,5L))
    c(solve(G,p0),l,r,mid,half)
  },numeric(10L)))
}

.np_udens_bb_guide_plan <- function(grid,h) {
  list(grid=as.double(grid),h=as.double(h),operator=.np_udens_bb_operator(grid,h),
       rules=list(`32`=.np_udens_q6_rule(32L),`48`=.np_udens_q6_rule(48L)))
}

.np_udens_bb_guide <- function(plan,ab,budget,order) {
  if (!is.numeric(ab) || length(ab)!=2L || any(!is.finite(ab)) || any(ab<=0) ||
      !isTRUE(budget %in% c(80,100)) || !isTRUE(order %in% c(32L,48L)))
    stop("invalid beta-binomial guide rule",call.=FALSE)
  grid <- plan$grid;h <- plan$h
  low <- -(budget+abs(lbeta(ab[1],ab[2])))/ab[1]
  high <- (budget+abs(lbeta(ab[1],ab[2])))/ab[2]
  physical <- sort(unique(c(grid,grid-sqrt(5)*h,grid+sqrt(5)*h)))
  physical <- physical[physical>0&physical<1]
  base <- sort(unique(c(low,0,qlogis(physical),high)))
  if (any(!is.finite(base)) || base[1]!=low || tail(base,1)!=high)
    stop("unresolved beta-binomial guide support",call.=FALSE)
  count <- ceiling(diff(base)/5)
  if (sum(count)>=4000L) stop("beta-binomial guide node budget exceeded",call.=FALSE)
  cuts <- sort(unique(unlist(lapply(seq_along(count),function(k)
    seq(base[k],base[k+1L],length.out=count[k]+1L)))))
  rule <- plan$rules[[as.character(order)]]
  half <- diff(cuts)/2;mid <- head(cuts,-1)+half
  t <- as.vector(outer(rule$x,half)+rep(mid,each=order))
  w <- as.vector(outer(rule$w,half))
  softplus <- function(t) pmax(t,0)+log1p(exp(-abs(t)))
  lv <- -softplus(-t);lq <- -softplus(t)
  prob <- exp(ab[1]*lv+ab[2]*lq-lbeta(ab[1],ab[2]))*w
  value <- .Call(C_np_plot_udens_bb_guide,grid,h,lv,lq,prob,plan$operator)
  if (any(!is.finite(value))) stop("nonfinite beta-binomial guide integral",call.=FALSE)
  list(value=value,mass=sum(prob),nodes=length(t),segments=length(half))
}

.np_udens_bb_factor <- function(q,T,AX,HI,m,n) {
  ans <- .Call(C_np_plot_udens_bb_factor,q,T,AX,HI,as.double(m),as.double(n))
  list(factor=attr(ans,"factor",exact=TRUE),mean=ans[,1L],variance=ans[,2L],
       correction=ans[,3L],moments=ans[,4:6,drop=FALSE])
}

# Reuse the unchanged complete analytic U5 formula, then represent its influence
# with certified positive weights. The expanded covariance remains an oracle.
.np_udens_bb_covariance <- function(z,x,grid,order) {
  .np_udens_u5_budget(length(x),length(grid),block=length(grid))
  j <- 0:z$k;A <- z$a+j;BB <- z$b+z$k-j
  weights <- rowMeans(vapply(x,function(t) dbinom(j,z$k,t),numeric(length(j))))
  moments <- z[["bb.log.moments",exact=TRUE]]
  if (is.null(moments)) moments <- .np_udens_bb_log_moments(c(z$a,z$b),z$k,weights)
  old <- .np_udens_u5_covariance(z,x,"epanechnikov",grid,order,4L,log.moments=moments)
  component <- exp(sweep(outer(old$lv,A)+outer(old$l1v,BB),2L,lbeta(A,BB),"-"))
  q <- drop(component %*% weights)*old$w
  m <- attr(moments,"mean",exact=TRUE)
  factor <- .np_udens_bb_factor(q,old$T,old$AX,old$HI,m,length(x))
  covariance <- tcrossprod(factor$factor)
  checks <- c(expanded=.np_udens_bb_check(covariance,old$expanded,"complete expanded covariance"),
              diagonal=.np_udens_bb_check(diag(covariance),factor$variance,"factor diagonal"))
  old$C <- covariance;old$psi_mean <- factor$mean;old$positive.factor <- factor$factor
  old$positive.checks <- list(refinement=checks,correction=factor$correction,moments=factor$moments)
  old
}

# Component-level preparation shared by guide matching and both representations
# of the complete BB covariance. Keep the small sum in paired precision until
# the final rounding; rounded component means alone are not the sum oracle.
.np_udens_bb_log_moments <- function(ab,k,weights) {
  if (!is.numeric(ab) || !is.null(dim(ab)) || !is.numeric(k) ||
      !is.null(dim(k)) || !is.numeric(weights) || !is.null(dim(weights)))
    stop("beta-binomial log moments require numeric vectors",call.=FALSE)
  .Call(C_np_plot_udens_bb_log_moments,as.double(ab),as.double(k),as.double(weights))
}
