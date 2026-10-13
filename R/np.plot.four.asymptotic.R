# New inference dispatch only. Existing fit/plot point owners are authoritative.
.np_plot_four_policy <- function(family, bws, gradient = NULL,
                                  view = "sample", proper = FALSE) {
  spec <- npConditionalRegEngineSpec(bws, where = "four-family target policy")
  positive <- identical(spec$reg.engine,"lp") && any(spec$degree.engine>0L)
  continuous <- !is.null(gradient) && bws[["ixcon",exact=TRUE]][gradient]
  route <- if(is.null(gradient))"level" else if(continuous)"continuous" else "categorical"
  if(identical(family,"qreg")) {
    target <- switch(route,level="inverse_cdf",continuous="leading_cdf_derivative",
      categorical=if(positive)"own_root_contrast" else "same_root_full_delta")
    return(list(route=route,target=target,proper=FALSE,
      scale=if(identical(target,"same_root_full_delta"))"joint_factor" else "incumbent"))
  }
  if(!identical(family,"conmode"))stop("unknown four-family target policy")
  normalize <- isTRUE(proper) && (route=="level" ||
    (route=="categorical" && (identical(view,"fixed") || positive)))
  list(route=route,target=if(route=="continuous")"raw_density_derivative" else
    if(route=="categorical")if(normalize)"proper_endpoint_contrast" else "raw_endpoint_contrast" else
      if(normalize)"proper_probability" else "raw_probability",
    proper=normalize,positive=positive,owner=if(identical(view,"fixed"))"conditional_hat" else "native_conditional",
    scale=if(normalize||route!="level")"joint_factor" else "incumbent")
}

.np_plot_four_factor_se <- function(A) {
  if(!is.matrix(A)||any(!is.finite(A)))stop("invalid joint influence factor")
  vapply(seq_len(nrow(A)),function(i) {
    if(all(A[i,]==0))return(0)
    scale<-max(abs(A[i,]));out<-scale*sqrt(sum((A[i,]/scale)^2))
    if(!is.finite(out)||out<=0).np_plot_joint_unavailable("joint uncertainty is not representable")
    out
  },numeric(1))
}

# Four-family Gaussian plot adapters. Existing estimators, marginal uncertainty,
# regression influence owners and Gaussian calibrators remain unchanged.
.np_plot_four_calibrate <- function(A, se, alpha, draws, target,
                                      exact.diagonal = FALSE) {
  if (!is.matrix(A) || any(!is.finite(A)))
    stop("invalid four-family influence factor", call. = FALSE)
  if (length(se) != nrow(A) || any(!is.finite(se)) || any(se < 0))
    .np_plot_joint_unavailable("marginal uncertainty is unavailable on the declared grid")
  norm <- numeric(nrow(A))
  zero <- logical(nrow(A))
  correlation.factor <- A
  for (i in seq_len(nrow(A))) {
    # A squared norm can underflow even when the influence is nonzero. Only
    # exact zero entries certify a degenerate row. Standardize before squaring.
    zero[i] <- all(A[i, ] == 0)
    if (zero[i]) next
    scale <- max(abs(A[i, ]))
    row <- A[i, ] / scale
    size <- sqrt(sum(row^2))
    norm[i] <- scale * size
    if (!is.finite(norm[i]) || norm[i] <= 0)
      stop("influence variance is not representable", call. = FALSE)
    correlation.factor[i, ] <- row / size
  }
  if (any(zero != (se == 0)))
    stop("influence and marginal uncertainty disagree about zero variance", call. = FALSE)
  if (exact.diagonal && any(abs(norm[!zero] / se[!zero] - 1) >
      256 * .Machine$double.eps))
    stop("joint influence does not reproduce the marginal variance", call. = FALSE)
  # This unit-row factor represents correlation only. Displayed widths retain
  # the separately justified marginal SE approximation. The unscaled factor
  # and its independent covariance oracle are never changed to fit a diagonal.
  cal <- if (nrow(A) <= 512L)
    .np_plot_gaussian_calibrate(tcrossprod(correlation.factor), alpha = alpha, draws = draws) else
    .np_plot_gaussian_factor_calibrate(correlation.factor, alpha = alpha, draws = draws)
  cal$target <- target
  cal$marginal.scale <- "incumbent asymptotic standard errors"
  cal$dependence <- "observation influence correlation"
  cal
}

.np_plot_four_interval <- function(interval, calibration, se, band.type) {
  if (is.null(calibration)) return(interval)
  joint.se <- calibration$simultaneous.se
  if (is.null(joint.se)) joint.se <- se
  simultaneous <- cbind(joint.se, joint.se, deparse.level = 0L) * calibration$critical
  if (identical(band.type, "all"))
    interval$all.err$simultaneous <- simultaneous else
    interval$err[, 1:2] <- simultaneous
  interval$calibration <- calibration
  interval
}

.np_plot_lsq_joint_local <- function(bws, x, ex, tau, se, gradient = NULL,
                                gradient.order = 1L, alpha = .05,
                                gaussian.draws = 9999L, estimate = NULL) {
  index <- .nplsqreg_match_plot_tau(tau, bws[["tau", exact = TRUE]])
  one <- if (length(bws[["tau", exact = TRUE]]) == 1L) bws else
    bws[["tau.bws", exact = TRUE]][[index]]
  if (is.null(one)) stop("missing per-tau location-scale state", call. = FALSE)
  reg <- one[["reg.bws", exact = TRUE]]
  deriv <- categorical <- NULL
  if (!is.null(gradient)) {
    continuous <- which(reg[["icon", exact = TRUE]])
    if (gradient %in% continuous) {
      deriv <- integer(length(continuous))
      deriv[match(gradient, continuous)] <- rep_len(gradient.order,
        length(continuous))[match(gradient, continuous)]
    } else categorical <- gradient
  }
  if (is.null(deriv)) {
    joint <- .np_plot_lsq_level_factor(reg, x, one[["qdat", exact = TRUE]],
      ex, categorical.index = categorical)
    if (is.null(estimate) || any(.np_plot_plreg_check_reconstruction(
        estimate, joint$estimate, se)))
      .np_plot_joint_unavailable("location-scale influence does not reconstruct the plotted target")
    influence <- joint$A
  } else {
    joint <- .np_plot_lsq_derivative_factor(reg,x,one[["qdat",exact=TRUE]],
      ex,coordinate=gradient,derivative.order=max(deriv),estimate=estimate,se=se)
    influence <- joint$A
  }
  .np_plot_four_calibrate(influence, se, alpha, gaussian.draws,
    "location-scale transformed-response regression; frozen pilot, delta and bandwidths")
}

.np_plot_quantile_interval <- function(interval, bws, x, y, ex, tau,
                                       estimate, se, gradient = NULL,
                                       gradient.order = 1L, alpha = .05,
                                       band.type = "all", gaussian.draws = 9999L) {
  calibration <- .np_plot_joint_try(if (inherits(bws, "lsqregressionbandwidth"))
    .np_plot_lsq_joint(bws, x, ex, tau, se, gradient, gradient.order,
      alpha, gaussian.draws, estimate = estimate) else
    .np_plot_qreg_joint(bws, x, y, ex, tau, estimate, se, gradient,
      alpha, gaussian.draws))
  .np_plot_four_interval(interval, calibration, se, band.type)
}

.np_plot_quantile_band_fields <- function(object, interval, gradient = NULL,
                                          tau.index = NULL) {
  if (is.null(tau.index))
    return(.np_plot_joint_band_fields(object, interval$calibration,
      interval$all.err, gradient))
  key <- if (is.null(gradient)) "quantile" else paste0("gradient", gradient)
  if (!is.null(interval$calibration))
    object$band.calibration[[key]][[tau.index]] <- interval$calibration
  if (!is.null(interval$all.err))
    object$band.errors[[key]][[tau.index]] <- interval$all.err
  object
}

.np_plot_four_budget <- function(n, m, p, width = 1L,
                                  max.bytes = 128 * 1024^2) {
  if (n < 2L || m < 1L || !is.finite(width) || width < 1L)
    stop("invalid four-family influence dimensions", call. = FALSE)
  tile <- min(m, max(1L, floor(4 * 1024^2 / (8 * as.double(n) * (12 + p)))))
  forecast <- 32 * as.double(n) * m + 64 * 1024^2 +
    8 * as.double(n) * tile * (18 + p) +
    64 * (as.double(n) + tile * (p + 2)) * width + 256 * width^2
  if (!is.finite(forecast) || forecast > max.bytes)
    .np_plot_joint_unavailable("four-family influence exceeds the plot memory budget")
  list(tile = tile, forecast.bytes = forecast)
}

.np_plot_four_conditional_admit <- function(bws, family) {
  spec <- npConditionalRegEngineSpec(bws, where = paste(family, "joint plot"))
  if (identical(spec$reg.engine, "lc") || all(spec$degree.engine == 0L)) 1L else
    npLpBasisNcol(spec$basis.engine, spec$degree.engine)
}

.np_plot_qreg_influence <- function(bws, x, y, ex, tau, quantile, se, gradient = NULL) {
  width <- .np_plot_four_conditional_admit(bws, "quantile")
  n <- nrow(x); m <- nrow(ex)
  budget <- .np_plot_four_budget(n, m, ncol(x), width)
  if (length(quantile) != m || any(!is.finite(quantile)) ||
      length(se) != m || any(!is.finite(se)) || any(se < 0))
    .np_plot_joint_unavailable("quantiles or marginal uncertainty are unavailable on the declared grid")
  if (ncol(y) != 1L || !is.numeric(y[[1L]]))
    .np_plot_joint_unavailable("regular quantile joint bands require one continuous response")
  # The inverter clamps to observed support when the requested probability
  # lies outside its bracket. A small residual alone does not make that root
  # regular. Metadata can be lost during plot row selection, so retain the
  # conservative endpoint certificate as well; never alter the fitted root.
  support <- range(y[[1L]], finite = TRUE)
  if (any(.npqreg_quantile_clamp(quantile) != "none") ||
      any(quantile <= support[1L] | quantile >= support[2L]))
    .np_plot_joint_unavailable("quantile joint bands require interior, unclamped roots")
  A <- matrix(NA_real_, m, n)
  density <- root <- reconstructed <- numeric(m)
  progress <- .np_plot_bootstrap_progress_begin(m, "Computing quantile plot influence", unit = "row")
  on.exit(.np_plot_progress_end(progress), add = TRUE)
  for (first in seq.int(1L, m, by = budget$tile)) {
    rows <- seq.int(first, min(m, first + budget$tile - 1L))
    ey <- setNames(data.frame(quantile[rows]), names(y))
    value <- .np_plot_four_conditional_factor(bws, x, y, ex[rows, , drop = FALSE], ey, cdf = TRUE, gradient = gradient)
    dens <- .np_plot_conditional_eval(bws = bws, xdat = x, ydat = y,
      exdat = ex[rows, , drop = FALSE], eydat = ey, cdf = FALSE,
      gradients = FALSE, se = FALSE, allow.external = TRUE)[["condens", exact = TRUE]]
    if (any(!is.finite(dens)) || any(dens <= 0))
      .np_plot_joint_unavailable("quantile joint bands require finite positive density slopes")
    level <- if (is.null(gradient)) value$estimate else value$level
    if (is.null(level)) level <- .np_plot_conditional_eval(bws, x, y,
      ex[rows, , drop = FALSE], ey, cdf = TRUE, gradients = FALSE,
      se = FALSE, allow.external = TRUE)[["condist", exact = TRUE]]
    residual <- level - tau
    # Inversion error must be negligible on the displayed uncertainty scale.
    # Do not refit, floor a slope, or accept a clamped tail root as regular.
    tolerance <- .01 * se[rows] * dens +
      max(128, 8 * as.double(n)) * .Machine$double.eps
    if (any(abs(residual) > tolerance))
      .np_plot_joint_unavailable("quantile inversion residual exceeds 1% of its asymptotic uncertainty")
    A[rows, ] <- -value$A / dens
    density[rows] <- dens; root[rows] <- residual
    reconstructed[rows] <- if (is.null(gradient)) quantile[rows] else -value$estimate / dens
    progress <- .np_plot_progress_tick(progress, max(rows))
  }
  list(A = A, density = density, root.residual = root, reconstructed = reconstructed,
    forecast.bytes = budget$forecast.bytes)
}

.np_plot_qreg_joint_local <- function(bws, x, y, ex, tau, estimate, se,
                                 gradient = NULL, alpha = .05,
                                 gaussian.draws = 9999L) {
  policy <- .np_plot_four_policy("qreg", bws, gradient)
  if (is.null(gradient)) {
    value <- .np_plot_qreg_influence(bws, x, y, ex, tau, estimate, se)
    target <- "inverse conditional CDF quantiles; selected smoothing, within one tau curve"
  } else {
    # The unchanged inverter owns each endpoint's quantile and clamp metadata.
    # Derivative rows retain the incumbent leading delta convention, requiring
    # shrinking regular neighborhoods and h_x^2/h_y -> 0 for its coverage law.
    quantile.fit <- function(z) npqreg(bws = bws, txdat = x, tydat = y,
      exdat = z, tau = tau, gradients = FALSE, se = TRUE)
    if (bws[["ixcon", exact = TRUE]][gradient]) {
      fit <- quantile.fit(ex)
      value <- .np_plot_qreg_influence(bws, x, y, ex, tau, fitted(fit),
        se(fit), gradient = gradient)
      target <- "conditional quantile continuous derivative; incumbent leading CDF-derivative delta approximation"
    } else if (identical(policy$target,"same_root_full_delta")) {
      fit <- quantile.fit(ex);q <- fitted(fit);root.se <- se(fit)
      support <- range(y[[1L]],finite=TRUE)
      if(any(.npqreg_quantile_clamp(q)!="none") ||
         any(!is.finite(c(q,root.se))) || any(root.se<0) ||
         any(q<=support[1L] | q>=support[2L]))
        .np_plot_joint_unavailable("quantile categorical full delta requires interior unclamped roots with finite uncertainty")
      value <- .np_plot_qreg_lc_same_root_delta(bws,x,y,ex,q,gradient,
        max.bytes=64*1024^2)
      if(any(abs(value$root-tau) > .01*root.se*value$primitives$Fy +
           max(128,8*as.double(nrow(x)))*.Machine$double.eps))
        .np_plot_joint_unavailable("quantile categorical full delta requires regular accurately inverted roots")
      value$reconstructed<-value$estimate
      value$root.residual<-value$root-tau
      value$primitives<-NULL
      target<-value$target
    } else {
      pair <- npCategoricalFirstDifferenceFrames(ex, gradient,
        where = "quantile joint plot")
      upper <- quantile.fit(pair$upper); lower <- quantile.fit(pair$lower)
      value <- .np_plot_qreg_influence(bws, x, y, pair$upper, tau,
        fitted(upper), se(upper))
      low <- .np_plot_qreg_influence(bws, x, y, pair$lower, tau,
        fitted(lower), se(lower))
      value$A <- value$A - low$A
      value$reconstructed <- value$reconstructed - low$reconstructed
      value$root.residual <- pmax(abs(value$root.residual), abs(low$root.residual))
      target <- "conditional quantile categorical contrast; paired own-root inverse CDF influence"
    }
    if (any(.np_plot_plreg_check_reconstruction(estimate, value$reconstructed, se)))
      .np_plot_joint_unavailable("quantile derivative influence does not reconstruct the plotted target")
  }
  joint.se <- if(identical(policy$scale,"joint_factor"))value$se else se
  cal <- .np_plot_four_calibrate(value$A, joint.se, alpha, gaussian.draws,
    target)
  if(identical(policy$scale,"joint_factor")) {
    cal$simultaneous.se<-joint.se
    cal$marginal.scale<-"new full-delta simultaneous uncertainty; existing pointwise and Bonferroni scales retained"
  }
  cal$maximum.root.residual <- max(abs(value$root.residual))
  cal$forecast.bytes <- value$forecast.bytes
  cal
}

.np_plot_conmode_raw_factor <- function(bws, x, y, ex, ey,
                                        gradient = NULL, view = "sample") {
  if (!identical(view, "fixed"))
    return(.np_plot_four_conditional_factor(bws, x, y, ex, ey,
      gradient = gradient))
  continuous <- which(bws[["ixcon", exact = TRUE]])
  s <- NULL
  if (!is.null(gradient)) {
    s <- integer(length(continuous)); s[match(gradient, continuous)] <- 1L
  }
  H <- .npcdhat_make_xhat_matrix(bws, x, ex, s = s)
  yb <- .npcdhat_make_ybw(bws, y)
  Y <- .npcdhat_make_kernel_matrix(yb, y, ey, operator = "normal")
  K <- .np_plot_density_matrix(.npcdhat_make_xkbw(bws, x), x, ex,
    rep("normal", ncol(x)))
  den <- rowSums(K)
  if (any(!is.finite(c(H, Y, K, den))) || any(den == 0))
    .np_plot_joint_unavailable("fixed-view class response map is unavailable")
  mu <- rowSums(K * Y) / den
  list(A = H * (Y - mu), estimate = rowSums(H * Y),
    target = "fixed-view conditional hat response map; local response moments")
}

.np_plot_conmode_delta_row <- function(raw, chosen, proper, tol) {
  K <- length(raw); jac <- numeric(K); jac[chosen] <- 1
  if (!proper) return(list(point = raw[chosen], jacobian = jac, zeros = integer()))
  if (any(!is.finite(raw))) .np_plot_joint_unavailable("nonfinite class probabilities")
  gamma <- function(k) { u <- k * .Machine$double.eps / 2; u / (1-u) }
  total <- sum(raw); margin <- max(-tol-min(raw), abs(total-1)-tol)
  roundoff <- if (is.infinite(tol)) 0 else gamma(K+4)*(sum(abs(raw))+1+abs(tol))
  if (abs(margin) <= roundoff)
    .np_plot_joint_unavailable("class probability normalizer is on a branch boundary")
  if (margin > 0) {
    point <- .npConmodeProjectSimplexVector(raw); active <- point > 0
    theta <- (sum(raw[active])-1)/sum(active)
    distance <- raw-theta
    bound <- gamma(sum(active)+3)*(sum(abs(raw[active]))+1)/sum(active) +
      gamma(2)*(abs(raw)+abs(theta))
    if (any(abs(distance) <= bound))
      .np_plot_joint_unavailable("class probability projection is on a face boundary")
    jac[] <- 0
    if (active[chosen]) {jac[active] <- -1/sum(active);jac[chosen] <- jac[chosen]+1}
    return(list(point=point[chosen],jacobian=jac,zeros=integer()))
  }
  active <- raw > 0; positive <- pmax(raw,0); total <- sum(positive)
  if (!is.finite(total) || total <= 0)
    .np_plot_joint_unavailable("class probability normalization has no positive mass")
  point <- positive[chosen]/total
  jac <- -point*as.numeric(active)/total
  if (active[chosen]) jac[chosen] <- jac[chosen]+1/total
  list(point=point,jacobian=jac,zeros=which(raw==0))
}

.np_plot_conmode_influence <- function(object, x, y, ex, level, levels,
                                        estimate = NULL, se = NULL,
                                        gradient = NULL, view = "sample") {
  bws <- object[["bws",exact=TRUE]];n<-nrow(x);m<-nrow(ex)
  width <- .np_plot_four_conditional_admit(bws,"class-probability")
  budget <- .np_plot_four_budget(n,m,ncol(x),width)
  if (8*as.double(m)*length(levels)*3+budget$forecast.bytes >128*1024^2)
    .np_plot_joint_unavailable("class support exceeds the plot memory budget")
  yy<-toFrame(y);names(yy)<-bws[["ynames",exact=TRUE]]
  chosen<-match(level,levels);if(is.na(chosen))stop("unknown class probability level")
  policy <- .np_plot_four_policy("conmode",bws,gradient,view,
    isTRUE(object[["proper.requested",exact=TRUE]]))
  if(!is.null(gradient)&&!bws[["ixcon",exact=TRUE]][gradient] &&
      (identical(view,"fixed") || policy$positive)) {
    endpoint.object <- object
    endpoint.object$proper.requested <- policy$proper
    pair<-npCategoricalFirstDifferenceFrames(ex,index=gradient,where="class joint plot")
    upper<-.np_plot_conmode_influence(endpoint.object,x,y,pair$upper,level,levels,view=view)
    lower<-.np_plot_conmode_influence(endpoint.object,x,y,pair$lower,level,levels,view=view)
    value<-list(A=upper$A-lower$A,estimate=upper$estimate-lower$estimate,
      forecast.bytes=budget$forecast.bytes,target="paired proper class endpoint contrast")
  } else {
    A<-matrix(NA_real_,m,n);point<-numeric(m)
    proper<-policy$proper
    tol<-object[["proper.info",exact=TRUE]][["tol",exact=TRUE]]
    if(is.null(tol))tol<-sqrt(.Machine$double.eps)
    progress<-.np_plot_bootstrap_progress_begin(m,"Computing class-probability plot influence",unit="row")
    on.exit(.np_plot_progress_end(progress),add=TRUE)
    for(first in seq.int(1L,m,by=budget$tile)) {
      rows<-seq.int(first,min(m,first+budget$tile-1L));z<-ex[rows,,drop=FALSE]
      make.y<-function(k)toFrame(.np_plot_conmode_level_factor(y,levels[k],length(rows)))
      if(!proper) {
        raw<-.np_plot_conmode_raw_factor(bws,x,yy,z,make.y(chosen),gradient,view)
        A[rows,]<-raw$A;point[rows]<-raw$estimate
      } else {
        raw<-matrix(NA_real_,length(rows),length(levels))
        for(k in seq_along(levels)) {
          ey<-make.y(k)
          raw[,k]<-if(identical(view,"fixed")) as.vector(npcdenshat(bws,txdat=x,
            tydat=yy,exdat=z,eydat=ey,y=rep.int(1,n),output="apply")) else
            .np_plot_conditional_eval(bws,x,yy,z,ey,cdf=FALSE,gradients=FALSE,
              se=FALSE,allow.external=TRUE)[["condens",exact=TRUE]]
        }
        maps<-lapply(seq_len(nrow(raw)),function(i).np_plot_conmode_delta_row(raw[i,],chosen,TRUE,tol))
        jac<-do.call(rbind,lapply(maps,`[[`,"jacobian"));tile<-matrix(0,length(rows),n)
        for(k in seq_along(levels)) {
          value<-.np_plot_conmode_raw_factor(bws,x,yy,z,make.y(k),view=view)
          for(i in seq_along(rows))if(k %in% maps[[i]]$zeros && any(value$A[i,]!=0))
            .np_plot_joint_unavailable("zero class probability lacks an exact zero-influence certificate")
          tile<-tile+value$A*jac[,k]
        }
        A[rows,]<-tile;point[rows]<-vapply(maps,`[[`,numeric(1),"point")
      }
      progress<-.np_plot_progress_tick(progress,max(rows))
    }
    value<-list(A=A,estimate=point,forecast.bytes=budget$forecast.bytes,
      target=if(is.null(gradient))"proper class probability on its stable local branch" else
        "raw conditional-density continuous derivative; incumbent conmode semantics")
  }
  value$marginal<-.np_plot_four_factor_se(value$A)
  value$policy<-policy
  value$target<-paste(policy$target,policy$owner,sep="; ")
  if(any(!is.finite(c(value$A,value$estimate,value$marginal))))
    .np_plot_joint_unavailable("class probability influence is nonfinite")
  if(!is.null(estimate)) {
    scale<-if(is.null(se))value$marginal else ifelse(is.finite(se),se,value$marginal)
    if(any(.np_plot_plreg_check_reconstruction(estimate,value$estimate,scale)))
      .np_plot_joint_unavailable("class probability influence does not reconstruct the plotted target")
  }
  value
}

 .np_plot_conmode_joint_columns_local <- function(dat, object, x, y, ex, level,
                                            levels, se, alpha, band,
                                            gaussian.draws) {
  gradient <- if(isTRUE(dat$gradients[1L]))match(dat$variable[1L],names(x)) else NULL
  target <- if(is.null(gradient))dat$probability else dat$effect
  view <- if(identical(dat$view[1L],"fixed"))"fixed" else "sample"
  value <- NULL
  cal <- .np_plot_joint_try({
    value <- .np_plot_conmode_influence(object,x,y,ex,level,levels,
      target,se,gradient,view)
    joint.se <- if(identical(value$policy$scale,"joint_factor"))value$marginal else se
    result <- if(band %in% c("simultaneous","all"))
      .np_plot_four_calibrate(value$A,joint.se,alpha,gaussian.draws,value$target) else
      list(method="normal marginal calibration", draws=0L,
        critical=qnorm(1-alpha/(2*if(identical(band,"bonferroni"))nrow(dat) else 1)),
        target=value$target, grid.size=nrow(dat))
    result$simultaneous.se <- joint.se
    result$marginal.scale <- if(identical(value$policy$scale,"joint_factor"))
      "new simultaneous-request delta standard errors; existing marginal bands retain their own scale" else
      "incumbent asymptotic standard errors"
    result$forecast.bytes <- value$forecast.bytes
    result
  })
  if(is.null(cal))return(dat)
  if(is.null(dat$center)) {
    # Gradient intervals are new; all their marginal/joint widths use the same
    # justified plot-only factor. Existing probability marginal fields survive.
    interval <- .np_plot_interval_payload(target,cal$simultaneous.se,
      plot.errors.method="asymptotic",plot.errors.alpha=alpha,plot.errors.type=band,
      plot.errors.center="estimate")
    dat$center<-interval$center;dat$stderr<-cal$simultaneous.se
    dat$lower<-target-interval$err[,1L];dat$upper<-target+interval$err[,2L]
    if(!is.null(interval$all.err))for(nm in names(interval$all.err)) {
      dat[[paste0(nm,".lower")]]<-target-interval$all.err[[nm]][,1L]
      dat[[paste0(nm,".upper")]]<-target+interval$all.err[[nm]][,2L]
    }
  }
  lower<-dat$center-cal$critical*cal$simultaneous.se
  upper<-dat$center+cal$critical*cal$simultaneous.se
  if(identical(band,"all")) {
    dat$simultaneous.lower<-lower;dat$simultaneous.upper<-upper
  } else if(identical(band,"simultaneous")) {dat$lower<-lower;dat$upper<-upper}
  attr(dat,"band.calibration")<-cal
  dat
}

.np_plot_copula_influence <- function(object, data, ex, se) {
  bws <- object[["bws", exact = TRUE]]
  if (identical(bws[["type", exact = TRUE]], "adaptive_nn")) {
    # New simultaneous inference owns its scale; legacy marginal fields may
    # remain NA. Native point owners and old interval requests are untouched.
    value <- tryCatch({
      helper <- if(identical(object[["target",exact=TRUE]],"density"))
        .np_plot_ann_bounded_density_factor else .np_plot_ann_bounded_distribution_factor
      helper(bws,data,ex,owner=environment(.np_plot_copula_influence))
    }, error=function(e) {
      prefixes <- c("np_plot_ann_unavailable: ", "np_plot_beta_numerical_unavailable: ",
        "np_plot_cdf_numerical_unavailable: ", "np_plot_nn_unavailable: ")
      which <- which(vapply(prefixes,function(p)startsWith(conditionMessage(e),p),logical(1)))
      if(length(which)) .np_plot_joint_unavailable(substring(conditionMessage(e),nchar(prefixes[which[1L]])+1L))
      stop(e)
    })
    estimate <- object[["copula",exact=TRUE]]
    roundoff <- max(128,8*as.double(nrow(data)))*.Machine$double.eps *
      pmax(1,abs(value$estimate),abs(estimate))
    if(length(estimate)!=nrow(ex)||any(!is.finite(value$estimate))||
       any(abs(value$estimate-estimate)>.01*value$simultaneous.se+roundoff))
      stop("ANN copula influence does not reproduce the plotted target",call.=FALSE)
    return(value)
  }
  if (!identical(bws[["type", exact = TRUE]], "fixed") &&
      identical(object[["target", exact = TRUE]], "density"))
    return(.np_plot_copula_nn_density_factor(object, data, ex, se,
      owner = environment(.np_plot_copula_influence)))
  n <- nrow(data); m <- nrow(ex); p <- ncol(data)
  budget <- .np_plot_four_budget(n, m, p)
  density <- identical(object[["target", exact = TRUE]], "density")
  if (length(se) != m || any(!is.finite(se)) || any(se < 0))
    .np_plot_joint_unavailable("copula uncertainty is unavailable on the declared grid")
  A <- matrix(NA_real_, m, n); fitted <- numeric(m)
  if (density) {
    marginal <- lapply(seq_len(p), function(j)
      .npcopula_marginal_bw(bws, data, j, target = "density", kbandwidth = TRUE))
    train <- lapply(seq_len(p), function(j) .npcopula_marginal_data(bws, data, j))
    evaluate <- lapply(seq_len(p), function(j) .npcopula_marginal_eval_data(bws, ex, j))
  }
  progress <- .np_plot_bootstrap_progress_begin(m, "Computing copula plot influence", unit = "row")
  on.exit(.np_plot_progress_end(progress), add = TRUE)
  for (first in seq.int(1L, m, by = budget$tile)) {
    rows <- seq.int(first, min(m, first + budget$tile - 1L))
    if (density) {
      product <- matrix(1, length(rows), n); total <- matrix(0, length(rows), n)
      for (j in seq_len(p)) {
        W <- .np_plot_density_matrix(marginal[[j]], train[[j]],
          evaluate[[j]][rows, , drop = FALSE], "normal")
        mass <- rowMeans(W)
        if (any(!is.finite(mass)) || any(mass == 0))
          .np_plot_joint_unavailable("copula density has a zero or nonfinite marginal denominator")
        W <- W / mass
        product <- product * W; total <- total + W
      }
      fitted[rows] <- rowMeans(product)
      phi <- product - fitted[rows] * (total - p + 1)
      A[rows, ] <- (phi - rowMeans(phi)) / n
    } else {
      value <- .np_plot_density_unconditional(bws, data, ex[rows, , drop = FALSE], cdf = TRUE)
      A[rows, ] <- value$A; fitted[rows] <- value$estimate
    }
    progress <- .np_plot_progress_tick(progress, max(rows))
  }
  estimate <- object[["copula", exact = TRUE]]
  roundoff <- max(128, 8 * as.double(n)) * .Machine$double.eps *
    pmax(1, abs(fitted), abs(estimate))
  if (length(estimate) != m || any(!is.finite(fitted)) ||
      any(abs(fitted - estimate) > .01 * se + roundoff))
    stop("copula influence does not reproduce the plotted target", call. = FALSE)
  list(A = A, estimate = fitted, forecast.bytes = budget$forecast.bytes,
    exact.diagonal = density)
}

.np_plot_copula_interval_local <- function(interval, object, data, ex, se,
                                      alpha, band.type, gaussian.draws) {
  cal <- .np_plot_joint_try({
    value <- .np_plot_copula_influence(object, data, ex, se)
    joint.se <- value$simultaneous.se
    if(is.null(joint.se)) joint.se <- se
    result <- .np_plot_four_calibrate(value$A, joint.se, alpha, gaussian.draws,
      if (!is.null(value$target)) value$target else if (value$exact.diagonal)
        "copula density ratio; frozen marginal inverse coordinates and bandwidths" else
        "copula distribution; frozen marginal inverse coordinates and bandwidths",
      exact.diagonal = value$exact.diagonal)
    result$forecast.bytes <- value$forecast.bytes
    result$simultaneous.se <- joint.se
    if(!is.null(value$simultaneous.se))
      result$marginal.scale <- "new simultaneous influence scale; incumbent marginal intervals retained"
    result
  })
  .np_plot_four_interval(interval, cal, se, band.type)
}


.np_plot_copula_interval <- function(...)
  .npRmpi_with_local_regression(.np_plot_copula_interval_local(...))

.np_plot_conmode_all_panel <- function(dat, ..., legend = TRUE, band = "all") {
  args <- list(...)
  # An unordered sample panel is an incumbent distributional boxplot, not a
  # curve at fixed covariates. Preserve that renderer and its interpretation.
  if (is.factor(dat$x) && !is.ordered(dat$x) &&
      !identical(dat$view[1L], "fixed"))
    return(do.call(.np_plot_conmode_panel, c(list(dat = dat), args)))
  bands <- if(identical(band,"all")) c("bonferroni", "simultaneous", "pointwise") else band
  columns <- if(identical(band,"all")) unlist(lapply(bands, function(b) paste0(b, c(".lower", ".upper")))) else c("lower","upper")
  if (!all(columns %in% names(dat)))
    stop("missing conmode joint interval columns", call. = FALSE)
  annotation <- args$annotation
  args$annotation <- NULL
  if (is.null(args$plot.user.args$ylim) && is.null(args$line.user.args$ylim))
    args$plot.user.args$ylim <- range(c(if(isTRUE(dat$gradients[1L]))dat$effect else dat$probability, unlist(dat[columns])), finite = TRUE)
  # Let the unchanged legacy panel own the central curve and axes. Supply no
  # legacy interval so it cannot overdraw one band with an ambiguous style.
  central <- dat
  central$lower <- central$upper <- NULL
  do.call(.np_plot_conmode_panel, c(list(dat = central), args))
  x <- if (is.factor(dat$x)) as.numeric(dat$x) else dat$x
  ord <- order(x)
  colors <- .np_plot_all_band_colors()
  # The accepted normal-band alias uses the same palette role as pointwise.
  colors[["pmzsd"]] <- colors[["pointwise"]]
  for (band in bands) {
    lower <- dat[[if(length(bands)==1L)"lower" else paste0(band, ".lower")]][ord]
    upper <- dat[[if(length(bands)==1L)"upper" else paste0(band, ".upper")]][ord]
    if (is.factor(dat$x) && !is.ordered(dat$x))
      graphics::segments(x[ord], lower, x[ord], upper, col = colors[[band]],
                         lty = .np_plot_lty("solid"), lwd = .np_plot_lwd("band_all_1d")) else {
      graphics::lines(x[ord], lower, col = colors[[band]],
                      lty = .np_plot_lty("solid"), lwd = .np_plot_lwd("band_all_1d"))
      graphics::lines(x[ord], upper, col = colors[[band]],
                      lty = .np_plot_lty("solid"), lwd = .np_plot_lwd("band_all_1d"))
    }
  }
  if (!is.null(annotation)) {
    annotation$continuous <- !is.factor(dat$x)
    .np_plot_draw_variability_annotation(annotation)
  }
  if(length(bands)>1L) .np_plot_draw_all_band_legend(legend = legend,
    plot.errors.method = "asymptotic", drawn = c(pointwise = TRUE, simultaneous = TRUE, bonferroni = TRUE))
  invisible(NULL)
}


# New inference only: numerical certificates accompany accepted beta maps.
.np_plot_lsq_precision_propagate <- function(mode,args) {
  tryCatch(.Call("C_np_plot_lsq_precision_propagate",mode,args,PACKAGE="npRmpi"),
    error=function(e) {
      if(startsWith(conditionMessage(e),"np_plot_precision_unavailable: "))
        .np_plot_joint_unavailable(sub("^np_plot_precision_unavailable: ","",conditionMessage(e)))
      stop(e)
    })
}

.np_plot_lsq_certified_residuals <- function(bws,x,q,rows) {
  n <- nrow(x)
  train <- .np_plot_lsq_gate_level_operator(bws,x,q)
  interval <- NULL
  seen <- logical(n)
  checked <- function(ii) {
    h <- train(ii)
    if(!is.matrix(h)||!identical(dim(h),c(length(ii),n))||any(!is.finite(h)))
      stop("invalid LSQ accepted training map",call.=FALSE)
    error <- attr(h,"operator.error",exact=TRUE)
    if(!is.null(error)) {
      bounds <- .np_plot_lsq_precision_propagate("residual",list(h,error,as.double(q),as.integer(ii)))
      if(is.null(interval))interval <<- matrix(NA_real_,n,2L)
      interval[ii,] <<- bounds
      seen[ii] <<- TRUE
    }
    h
  }
  # This protected helper remains byte-for-byte unchanged. Its actual returned
  # residual, rather than an interval center, is used to form every factor.
  result <- .np_plot_response_residuals(bws,x,q,rows,training.operator=checked)
  if(any(!is.finite(result$residual)))
    .np_plot_joint_unavailable("LSQ normalized residuals are not finite")
  if(!is.null(interval) && (!all(seen)||any(!is.finite(interval))))
    stop("incomplete LSQ residual precision certificate",call.=FALSE)
  list(residual=result$residual,interval=interval)
}

.np_plot_lsq_map_difference <- function(upper,lower,same.endpoint=NULL) {
  eu <- attr(upper,"operator.error",exact=TRUE)
  el <- attr(lower,"operator.error",exact=TRUE)
  if(is.null(eu)&&is.null(el))return(upper-lower)
  if(is.null(eu)||is.null(el))stop("inconsistent LSQ endpoint precision certificates")
  if(is.null(same.endpoint))same.endpoint <- rep(FALSE,nrow(upper))
  if(!is.logical(same.endpoint)||length(same.endpoint)!=nrow(upper)||anyNA(same.endpoint))
    stop("invalid LSQ endpoint identity certificate")
  # An identical declared endpoint is the same random operator. Numerical
  # equality alone never establishes this certificate; the caller supplies
  # exact complete-frame identity and both accepted maps must also agree.
  if(any(upper[same.endpoint,,drop=FALSE]!=lower[same.endpoint,,drop=FALSE]))
    .np_plot_joint_unavailable("identical LSQ endpoints have inconsistent accepted maps")
  h <- error <- matrix(0,nrow(upper),ncol(upper))
  keep <- which(!same.endpoint)
  if(length(keep)) {
    out <- .np_plot_lsq_precision_propagate("difference",list(upper[keep,,drop=FALSE],
      eu[keep,,drop=FALSE],lower[keep,,drop=FALSE],el[keep,,drop=FALSE]))
    h[keep,] <- out$value
    error[keep,] <- out$error
  }
  attr(h,"operator.error") <- error
  h
}

.np_plot_lsq_check_factor_precision <- function(map,residual.interval,actual) {
  error <- attr(map,"operator.error",exact=TRUE)
  if(is.null(error)&&is.null(residual.interval))return(invisible(NULL))
  if(is.null(error)||is.null(residual.interval))
    stop("inconsistent LSQ factor precision certificates",call.=FALSE)
  proof <- .np_plot_lsq_precision_propagate("factor",list(map,error,residual.interval,actual))
  if(any(!proof$certified))
    .np_plot_joint_unavailable("accepted LSQ factor cannot certify requested floating-point precision")
  invisible(proof$relative.error)
}

.np_plot_lsq_level_factor <- function(bws, x, q, ex, categorical.index = NULL,
                                    block.bytes = 4 * 1024^2,
                                    max.bytes = 128 * 1024^2) {
  n <- nrow(x); m <- nrow(ex)
  if (!is.data.frame(x) || !is.data.frame(ex) || n < 2L || m < 1L ||
      !identical(names(x), names(ex)) || length(q) != n ||
      !is.numeric(q) || any(!is.finite(q)) || anyNA(x) || anyNA(ex))
    stop("invalid LSQ level-factor data", call. = FALSE)
  if (length(block.bytes) != 1L || length(max.bytes) != 1L ||
      !is.finite(block.bytes) || !is.finite(max.bytes) ||
      block.bytes <= 0 || max.bytes <= 0)
    stop("invalid LSQ level-factor budget", call. = FALSE)
  spec <- npValidatedConditionalRegSpec(bws,
    where = "LSQ actual level factor", ncon.field = "ncon")
  width <- if (identical(spec$regtype.engine, "lc")) 1 else
    npLpBasisNcol(spec$basis.engine, spec$degree.engine)
  # Two evaluation endpoint operators may coexist. Account for their compact
  # projections/bases, training operator, final factor, row temporaries and the
  # unchanged native owner's setup before allocating any estimator workspace.
  fixed <- 64 * 1024^2 + 64 * (as.double(n) + m) * (width + ncol(x) + 1) +
    272 * width^2 + 64 * as.double(n) * m + 64 * as.double(n)
  # The unchanged native ANN owner retains its own 64-row training tile;
  # requested R block size cannot make that private native allocation smaller.
  if (identical(bws[["type", exact = TRUE]], "adaptive_nn"))
    fixed <- fixed + 8 * as.double(n) * min(n, 64L)
  available <- min(block.bytes, max.bytes - fixed)
  rows <- min(n, m, floor(available / (160 * as.double(n))))
  if (!is.finite(rows) || rows < 1L)
    .np_plot_joint_unavailable("LSQ accepted-map factor exceeds the memory budget")
  if (!is.null(categorical.index) &&
      (length(categorical.index) != 1L || !is.finite(categorical.index) ||
       categorical.index != as.integer(categorical.index) ||
       categorical.index < 1L || categorical.index > ncol(x) ||
       !is.factor(x[[categorical.index]])))
    stop("invalid LSQ categorical factor coordinate", call. = FALSE)

  residual <- .np_plot_lsq_certified_residuals(bws,x,q,rows)
  z <- residual$residual
  if (is.null(categorical.index)) {
    upper <- .np_plot_lsq_gate_level_operator(bws, x, q, ex)
    lower <- NULL
    same.endpoint <- NULL
  } else {
    frames <- npCategoricalFirstDifferenceFrames(ex, categorical.index,
      "LSQ actual level factor")
    same.endpoint <- Reduce(`&`,Map(function(a,b)a==b,frames$upper,frames$lower))
    upper <- .np_plot_lsq_gate_level_operator(bws, x, q, frames$upper)
    lower <- .np_plot_lsq_gate_level_operator(bws, x, q, frames$lower)
  }
  A <- matrix(0, m, n); fitted <- numeric(m)
  for (first in seq.int(1L, m, by = rows)) {
    ii <- seq.int(first, min(m, first + rows - 1L))
    h <- upper(ii)
    if (!is.null(lower)) h <- .np_plot_lsq_map_difference(h,lower(ii),same.endpoint[ii])
    if (!is.matrix(h) || !identical(dim(h), c(length(ii), n)) ||
        any(!is.finite(h))) stop("invalid accepted evaluation map", call. = FALSE)
    fitted[ii] <- drop(h %*% q)
    tile <- sweep(h, 2L, z, "*")
    .np_plot_lsq_check_factor_precision(h,residual$interval,tile)
    A[ii, ] <- tile
  }
  list(A = A, residual = z, estimate = fitted,
       block.rows = rows, forecast.bytes = fixed + 160 * n * rows,
       target = "frozen LSQ transformed response; accepted native level maps")
}


.np_plot_copula_nn_density_factor <- function(object, data, ex, se, owner,
                                      max.bytes = 128 * 1024^2,
                                      tile.rows = NULL) {
  if (!is.environment(owner)) stop("owner must be the exact package/source environment")
  f <- function(name) get(name, envir = owner, mode = "function", inherits = TRUE)
  unavailable <- f(".np_plot_joint_unavailable")
  bws <- object[["bws", exact = TRUE]]
  type <- bws[["type", exact = TRUE]]
  if (!identical(object[["target", exact = TRUE]], "density") ||
      !type %in% c("generalized_nn", "adaptive_nn"))
    stop("NN copula factor requires an NN copula-density object")
  if (!is.data.frame(data) || !is.data.frame(ex) ||
      !identical(names(data), bws[["xnames", exact = TRUE]]) ||
      !identical(names(ex), names(data)) || anyNA(data) || anyNA(ex))
    stop("complete data/evaluation frames must follow the fitted variable order")
  n <- nrow(data); m <- nrow(ex); p <- ncol(data)
  pc <- sum(bws[["icon", exact = TRUE]])
  if (n < 2L || m < 1L || p < 1L || !is.finite(pc) || pc < 0 || pc > p)
    stop("invalid NN copula influence dimensions")
  if (!is.numeric(se) || length(se) != m || any(!is.finite(se)) || any(se < 0))
    unavailable("NN copula marginal uncertainty is unavailable on the declared grid")
  incumbent <- object[["copula", exact = TRUE]]
  if (!is.numeric(incumbent) || length(incumbent) != m || any(!is.finite(incumbent)))
    unavailable("NN copula density is unavailable on the declared grid")
  if (length(max.bytes) != 1L || !is.numeric(max.bytes) ||
      !is.finite(max.bytes) || max.bytes <= 0)
    stop("NN copula factor memory budget must be positive and finite")

  # A, correlation/copy allowance, unchanged Gaussian workspace, bounded ANN
  # normal/derivative/support-face tiles, and O(n*p) geometry. No class tensor.
  tile <- min(m, max(1, floor(4 * 1024^2 /
                             (8 * as.double(n) * (16 + 2 * pc)))))
  if (!is.null(tile.rows)) {
    if (!is.numeric(tile.rows) || length(tile.rows) != 1L ||
        !is.finite(tile.rows) || tile.rows != floor(tile.rows) ||
        tile.rows < 1 || tile.rows > tile)
      stop("tile.rows may only reduce the budgeted positive integer tile")
    tile <- tile.rows
  }
  forecast <- 24 * as.double(n) * m + 64 * 1024^2 +
    8 * as.double(n) * tile * (16 + 2 * pc) +
    256 * as.double(n) * (p + pc + 1) + 128 * as.double(m) * (p + 1)
  if (!is.finite(forecast) || forecast > max.bytes)
    unavailable("NN copula influence exceeds the plot memory budget")

  numerator.owner <- f(".np_plot_density_unconditional")
  marginal.bws <- lapply(seq_len(p), function(j)
    f(".npcopula_marginal_bw")(bws, data, j, target = "density", kbandwidth = TRUE))
  train <- lapply(seq_len(p), function(j) f(".npcopula_marginal_data")(bws, data, j))
  evaluate <- lapply(seq_len(p), function(j) f(".npcopula_marginal_eval_data")(bws, ex, j))
  A <- matrix(NA_real_, m, n)
  estimate <- denominator <- numeric(m)
  numerator.targets <- character()
  minimum.marginal <- Inf
  progress <- f(".np_plot_bootstrap_progress_begin")(
    m, "Computing nearest-neighbor copula plot influence", unit = "row")
  on.exit(f(".np_plot_progress_end")(progress), add = TRUE)
  for (first in seq.int(1L, m, by = tile)) {
    rows <- seq.int(first, min(m, first + tile - 1L))
    # Same primitives and multiplication order as the incumbent denominator,
    # but validate BEFORE its final NZD(). A floored zero is not positive proof.
    D <- rep.int(1.0, length(rows))
    for (j in seq_len(p)) {
      dj <- f(".np_ksum_unconditional_eval_exact")(
        xdat = train[[j]], exdat = evaluate[[j]][rows, , drop = FALSE],
        bws = marginal.bws[[j]], operator = "normal")
      if (length(dj) != length(rows) || any(!is.finite(dj)) || any(dj <= 0))
        unavailable("NN copula density requires strictly positive finite raw marginal densities")
      minimum.marginal <- min(minimum.marginal, dj)
      D <- D * dj
      if (any(!is.finite(D)) || any(D <= 0))
        unavailable("NN copula marginal-density product is not representable and positive")
    }
    value <- numerator.owner(bws, data, ex[rows, , drop = FALSE], cdf = FALSE)
    if (!is.matrix(value$A) || !identical(dim(value$A), c(length(rows), n)) ||
        any(!is.finite(value$A)) || length(value$estimate) != length(rows) ||
        any(!is.finite(value$estimate)))
      stop("unchanged NN numerator returned malformed/nonfinite influence")
    expected.target <- if (identical(type, "adaptive_nn") && pc > 0)
      "same-sample ANN random-radius influence" else
      "kernel empirical mass; fixed or frozen-query geometry"
    if (!identical(value$target, expected.target))
      stop("NN copula numerator owner target changed; review the binding")
    # R's column-major recycling divides every row by its own D. No recentering,
    # scaling to fit an SE, independent marginal draws or donor permutation.
    A[rows, ] <- value$A / D
    estimate[rows] <- value$estimate / D
    denominator[rows] <- D
    numerator.targets <- union(numerator.targets, value$target)
    value <- NULL
    progress <- f(".np_plot_progress_tick")(progress, max(rows))
  }
  if (any(!is.finite(A)) || any(!is.finite(estimate)))
    unavailable("NN copula factor division is not representable")
  # Preserve the maintained fit-reproduction predicate; do not manufacture a
  # diagonal equality between different incumbent asymptotic approximations.
  roundoff <- max(128, 8 * as.double(n)) * .Machine$double.eps *
    pmax(1, abs(estimate), abs(incumbent))
  if (any(abs(estimate - incumbent) > .01 * se + roundoff))
    stop("NN copula factor does not reproduce the plotted target within 1% of its SE")
  list(A = A, estimate = estimate, denominator = denominator,
       exact.diagonal = FALSE, numerator.target = numerator.targets,
       target = "NN copula density; frozen inverse coordinates and marginal-density denominator; incumbent numerator variability",
       marginal.scale = "incumbent NN numerator asymptotic SE divided by the same positive denominator",
       convention = list(inverse.coordinates = "frozen", denominator = "frozen",
                         numerator = numerator.targets, continuous.dimensions = pc,
                         leading.rate.interpretation = if (pc >= 2L)
                           "additional regular shrinking-NN assumptions and H/h_j -> 0 required" else
                           "working variability only; multivariate denominator-negligibility rate not claimed"),
       forecast.bytes = forecast, tile.rows = tile,
       order = list(evaluation = seq_len(m), observation = seq_len(n)),
       denominator.certificate = list(raw.positive = TRUE,
                                      minimum.marginal = minimum.marginal,
                                      minimum.product = min(denominator)))
}


.np_plot_four_conditional_factor <- function(bws, x, y, ex, ey, cdf = FALSE,
                                             gradient = NULL, order = 1L) {
  spec <- npConditionalRegEngineSpec(bws, where = "four-family conditional influence")
  if (identical(spec$reg.engine, "lc") || all(spec$degree.engine == 0L))
    return(.np_plot_density_conditional(bws, x, y, ex, ey, cdf,
      gradient = gradient, order = order))
  n <- nrow(x); m <- nrow(ex)
  if (n < 2L || m < 1L || nrow(y) != n || nrow(ey) != m)
    stop("unaligned conditional influence data", call. = FALSE)
  if (!is.null(gradient) && !bws[["ixcon", exact = TRUE]][gradient])
    stop("categorical influence requires paired complete endpoint fits", call. = FALSE)
  width <- npLpBasisNcol(spec$basis.engine, spec$degree.engine)
  value <- .np_conditional_eval_selected(bws, x, y, ex, ey, cdf = cdf,
    gradients = !is.null(gradient), categorical.effects = FALSE,
    gradient.order = order, se = TRUE, allow.external = TRUE,
    projection.request = as.integer(width))
  projection <- attr(value, ".np.conditional.projection", exact = TRUE)
  if (!is.list(projection) || !identical(dim(projection$responses), c(n, m)))
    stop("missing accepted conditional response/projection", call. = FALSE)
  A <- matrix(NA_real_, m, n); reconstructed <- center <- numeric(m)
  target <- if (is.null(gradient)) fitted(value) else gradients(value)[, gradient]
  marginal <- if (is.null(gradient)) se(value) else gradients(value, se = TRUE)[, gradient]
  # Request scope is one actual conditional row at a time. Its exact gate,
  # physical bandwidth, donor order and accepted adjoint are retained natively.
  for (j in seq_len(m)) {
    constant <- projection$constant[j] == 1L
    K <- if (constant) rep.int(1, n) else projection$weights[, j]
    if (any(!is.finite(K)) || !is.finite(sum(K)) || sum(K) == 0)
      .np_plot_joint_unavailable("conditional influence kernel mass is unavailable")
    if (is.null(gradient)) P <- projection$projection[j, , drop = FALSE] else {
      coordinate <- match(gradient, which(bws[["ixcon", exact = TRUE]]))
      if (constant) {
        H <- .Call("C_np_plot_alllarge_derivative", projection$basis,
          as.matrix(x[, bws[["ixcon", exact = TRUE]], drop = FALSE]),
          as.matrix(ex[j, bws[["ixcon", exact = TRUE]], drop = FALSE]),
          as.integer(spec$degree.engine), isTRUE(spec$bernstein.engine),
          as.integer(match(spec$basis.engine, c("additive", "glp", "tensor")) - 1L),
          as.integer(coordinate), as.integer(order), projection$donor.order,
          projection$projection[j, , drop = FALSE], PACKAGE = "npRmpi")$weights
      } else {
      slots <- which(projection$coordinate[, j] == coordinate &
        projection$order[, j] == order)
      if (length(slots) != 1L || slots > projection$captured[j])
        .np_plot_joint_unavailable("accepted conditional derivative direction is missing")
      col <- (j - 1L) * nrow(projection$coordinate) + slots
      P <- matrix(projection$adjoints[, col], nrow = 1L)
      }
    }
    if (is.null(gradient) || !constant) H <- .Call("C_np_plot_projection_apply", projection$basis, P,
      if (constant) NULL else matrix(K, ncol = 1L), PACKAGE = "npRmpi")
    Z <- projection$responses[, j]
    scale <- exp(projection$response.logscale[j])
    if (any(!is.finite(Z)) || !is.finite(scale))
      .np_plot_joint_unavailable("conditional response influence cannot be represented")
    mu <- if (all(Z == Z[1L])) Z[1L] else sum(K * Z) / sum(K)
    reconstructed[j] <- sum(H * Z) * scale
    center[j] <- mu * scale
    A[j, ] <- H * (Z - mu) * scale
  }
  if (any(!is.finite(A)) || any(!is.finite(target)) ||
      any(!is.finite(marginal)) || any(marginal < 0))
    .np_plot_joint_unavailable("accepted conditional uncertainty is unavailable")
  if (any(.np_plot_plreg_check_reconstruction(target, reconstructed, marginal)))
    .np_plot_joint_unavailable("conditional influence does not reconstruct the accepted target")
  list(A = A, estimate = target, reconstructed = reconstructed,
    level = fitted(value), center = center, marginal = marginal,
    target = "accepted conditional response map; local conditional response moments")
}

# New-request bridge for the new fixed native dispatcher. No protected function
# is replaced, traced, rebound, or edited.
.np_plot_lsq_gate_budget <- function(n, m, p, parent.rows, max.bytes, width=0) {
  cells <- as.double(n)*m
  forecast <- 64*cells + 64*(as.double(n)+parent.rows)*(p+width+1) +
    256*as.double(width)^2
  if(length(max.bytes)!=1L || !is.finite(max.bytes) || max.bytes<=0 ||
     !is.finite(cells) || cells>.Machine$integer.max ||
     !is.finite(forecast) || forecast>max.bytes)
    .np_plot_joint_unavailable('LSQ gate scope exceeds its memory or native dimension budget')
  invisible(forecast)
}

.np_plot_lsq_gate_geometry <- function(bws, x, ex = NULL, max.bytes=128*1024^2) {
  training <- is.null(ex)
  x <- toFrame(x)
  if (training) ex <- x else ex <- toFrame(ex)
  if (!(x %~% ex) || !nrow(x) || !nrow(ex) || anyNA(x) || anyNA(ex))
    stop('invalid LSQ gate geometry', call. = FALSE)
  .np_plot_lsq_gate_budget(nrow(x),0,ncol(x),nrow(ex),max.bytes)
  x <- adjustLevels(x, bws$xdati)
  ex <- adjustLevels(ex, bws$xdati, allowNewCells = TRUE)
  npKernelBoundsCheckEval(ex, bws$icon, bws$ckerlb, bws$ckerub,
                          argprefix = 'cker')
  split <- function(dat) {
    z <- toMatrix(dat)
    out <- lapply(list(uno=bws$iuno, ord=bws$iord, con=bws$icon),
      function(j) {
        a <- z[, j, drop=FALSE]
        storage.mode(a) <- 'double'
        a
      })
    out
  }
  bw <- .np_physical_bandwidth(bws)
  list(train=split(x), eval=split(ex),
       physical.bandwidth=as.double(c(bw[bws$icon],bw[bws$iuno],bw[bws$iord])),
       support=.np_native_categorical_support(bws),
       train.index=seq_len(nrow(x)), eval.index=seq_len(nrow(ex)))
}

.np_plot_lsq_gate_indices <- function(ii, n) {
  if (!is.integer(ii) || !length(ii) || anyNA(ii) || any(ii<1L | ii>n))
    stop('invalid LSQ gate row map', call.=FALSE)
  ii
}

# Main kernel rows use physical bandwidths and divided returned weights for
# every bandwidth type, matching the accepted projection reconstruction.
.np_plot_lsq_gate_weights <- function(bws, geometry, index=geometry$eval.index,
                             training=FALSE, package='npRmpi',max.bytes=128*1024^2) {
  n <- nrow(geometry$train$con); m <- nrow(geometry$eval$con)
  index <- .np_plot_lsq_gate_indices(index, m)
  .np_plot_lsq_gate_budget(n,length(index),length(bws$bw),m,max.bytes)
  kb <- if (isa(bws,'kbandwidth')) bws else kbandwidth(bws)
  opt <- list(n,length(index),kb$nuno,kb$nord,kb$ncon,SF_ARB,
    switch(kb$type,fixed=BW_FIXED,generalized_nn=BW_GEN_NN,adaptive_nn=BW_ADAP_NN),
    if(isTRUE(getOption('np.messages'))) IO_MIN_FALSE else IO_MIN_TRUE,
    switch(kb$ckertype,gaussian=CKER_GAUSS+kb$ckerorder/2-1,
      epanechnikov=CKER_EPAN+kb$ckerorder/2-1,uniform=CKER_UNI,beta=CKER_COORDINATE),
    switch(kb$ukertype,aitchisonaitken=UKER_AIT,liracine=UKER_LR),
    switch(kb$okertype,wangvanryzin=OKER_WANG,liracine=OKER_LR,
      nliracine=OKER_NLR,racineliyan=OKER_RLY),
    FALSE,FALSE,TRUE,attr(kb$xmcv,'num.row'),0L,0L,
    npDoTreeOrCategoricalCompress(ncon=kb$ncon,ncat=kb$nuno+kb$nord,bws=kb),
    TRUE,PERMUTATION_OPERATORS[['none']],FALSE,FALSE,
    suppress.parallel=TRUE)
  opt <- as.integer(c(opt,npContinuousKernelDescriptorOptions(kb),list(
    divide.returned.kernel.weights=TRUE,
    categorical.compress=identical(kb$ckertype,'beta') && kb$nuno+kb$nord>0L &&
      npStrictLogicalOption('np.categorical.compress',TRUE),
    tree.outer.blas=FALSE)))
  if (training && identical(kb$type,'generalized_nn') && kb$ncon>0L) {
    if (!identical(geometry$train,geometry$eval))
      stop('training kernel occurrence map requires identical geometry',call.=FALSE)
    attr(opt,'np.eval.train.index') <- index-1L
  }
  lohi <- npKernelBoundsMarshal(kb$ckerlb[kb$icon],kb$ckerub[kb$icon])
  t <- geometry$train
  e <- lapply(geometry$eval,function(a) a[index,,drop=FALSE])
  args <- list(t$uno,t$ord,t$con,as.double(0),as.double(0),
    e$uno,e$ord,e$con,geometry$physical.bandwidth,as.double(kb$xmcv),
    as.double(attr(kb$xmcv,'pad.num')),rep.int(as.integer(ALL_OPERATORS[['normal']]),
      kb$ncon+kb$nuno+kb$nord),opt,as.double(1),as.integer(length(index)),
    0L,as.integer(as.double(n)*length(index)),as.double(lohi$lb),as.double(lohi$ub))
  geometry$eval.index <- index
  out <- .np_with_nn_radius_context(
    .Call('C_np_plot_gate_dispatch','kernelsum',args,geometry,PACKAGE=package),
    continuous.names=bws[['xnames',exact=TRUE]][bws[['icon',exact=TRUE]]],
    bws=bws,ntrain=n)
  matrix(as.double(out[['kernel.weights']]),nrow=n,ncol=length(index))
}

# Whole-evaluation derivative maps preserve the native derivative owner and
# ordering. In particular, this does not use public hat higher-order dispatch.
.np_plot_lsq_gate_derivative <- function(bws,x,ex=NULL,s,package='npRmpi',max.bytes=128*1024^2) {
  geometry <- .np_plot_lsq_gate_geometry(bws,x,ex,max.bytes=max.bytes)
  t <- geometry$train; e <- geometry$eval
  spec <- npValidatedConditionalRegSpec(bws,where='LSQ scoped derivative',ncon.field='ncon')
  width <- if(npIsCanonicalLp0Spec(spec,bws$ncon)) 1 else
    npLpBasisNcol(spec$basis.engine,spec$degree.engine)
  .np_plot_lsq_gate_budget(nrow(t$con),nrow(e$con),length(bws$bw),nrow(e$con),max.bytes,width)
  if (!is.integer(s) || length(s)!=bws$ncon || anyNA(s) || any(s<0L))
    stop('invalid LSQ derivative vector',call.=FALSE)
  npCheckRegressionDesignCondition(
    reg.code=npRegtypeToC(regtype='lp',degree=spec$degree.engine,ncon=bws$ncon,
      context='LSQ scoped derivative')$code,
    xcon=t$con,basis=spec$basis.engine,degree=spec$degree.engine,
    bernstein.basis=spec$bernstein.basis.engine,where='LSQ scoped derivative')
  beta <- identical(bws[['ckertype',exact=TRUE]],'beta')
  lohi <- npKernelBoundsMarshal(bws$ckerlb[bws$icon],bws$ckerub[bws$icon])
  args <- list(t$uno,t$ord,t$con,e$uno,e$ord,e$con,
    matrix(double(),nrow(t$con),0L),geometry$physical.bandwidth,
    as.integer(switch(bws$type,fixed=BW_FIXED,generalized_nn=BW_GEN_NN,adaptive_nn=BW_ADAP_NN)),
    as.integer(if(beta) npContinuousKernelCode(bws) else switch(bws$ckertype,
      gaussian=CKER_GAUSS+bws$ckerorder/2-1,
      epanechnikov=CKER_EPAN+bws$ckerorder/2-1,uniform=CKER_UNI)),
    as.integer(switch(bws$ukertype,aitchisonaitken=UKER_AIT,liracine=UKER_LR)),
    as.integer(switch(bws$okertype,wangvanryzin=OKER_WANG,liracine=OKER_LR,racineliyan=OKER_RLY)),
    !beta && identical(.npreg_fit_tree_code(bws,ncon=bws$ncon,ncat=bws$nuno+bws$nord),DO_TREE_YES),
    as.integer(spec$degree.engine),if(any(s>0L)) s else integer(),
    as.integer(isTRUE(spec$bernstein.basis.engine)),as.integer(npLpBasisCode(spec$basis.engine)),
    as.integer(if(beta) CKER_FAMILY_BETA else CKER_FAMILY_LEGACY),as.integer(bws$ckerorder),
    as.double(lohi$lb),as.double(lohi$ub),
    as.integer(beta && bws$nuno+bws$nord>0L && npStrictLogicalOption('np.categorical.compress',TRUE)),
    TRUE,is.null(ex),FALSE,FALSE,.np_native_categorical_support(bws))
  on.exit(tryCatch(.Call('C_np_reset_native_estimator_state',PACKAGE=package),
    error=function(e) NULL),add=TRUE)
  .np_with_nn_radius_context(
    .Call('C_np_plot_gate_dispatch','lp_apply_ctx',args,geometry,PACKAGE=package),
    continuous.names=bws[['xnames',exact=TRUE]][bws[['icon',exact=TRUE]]],
    bws=bws,ntrain=nrow(t$con),leave.one.out=FALSE)
}

# Accepted level adjoints and a scoped reconstruction of their actual kernel
# weights. Each closure retains the full parent evaluation geometry.
.np_plot_lsq_gate_level_operator <- function(bws,x,q,ex=NULL,package='npRmpi',max.bytes=128*1024^2) {
  geometry <- .np_plot_lsq_gate_geometry(bws,x,ex,max.bytes=max.bytes)
  spec <- npValidatedConditionalRegSpec(bws,where='LSQ scoped level map',ncon.field='ncon')
  scalar <- npIsCanonicalLp0Spec(spec,bws$ncon) || bws$ncon<1L
  width <- if(scalar) 1 else npLpBasisNcol(spec$basis.engine,spec$degree.engine)
  .np_plot_lsq_gate_budget(nrow(x),0,ncol(x),nrow(geometry$eval$con),max.bytes,width)
  if (scalar) return(function(index) {
    kw <- .np_plot_lsq_gate_weights(bws,geometry,index,training=is.null(ex),package=package,max.bytes=max.bytes)
    .np_lc_hat_normalize(kw,colSums(kw))
  })
  fit <- .np_regression_direct(bws,txdat=x,tydat=q,exdat=ex,fit.projection=TRUE)
  ingredients <- fit$fit.projection
  if (!is.list(ingredients) || length(ingredients)!=3L)
    stop('missing LSQ accepted projection',call.=FALSE)
  function(index) {
    index <- .np_plot_lsq_gate_indices(index,nrow(geometry$eval$con))
    .np_plot_lsq_gate_budget(nrow(geometry$train$con),length(index),length(bws$bw),
      nrow(geometry$eval$con),max.bytes,width)
    if(identical(bws[['ckertype',exact=TRUE]],'beta')) {
      certificate <- .np_plot_lsq_beta_moment_weights(bws,x,q,ingredients[[1L]],ex,index,
        max.bytes=max.bytes,package=package,accepted.projection=ingredients[[2L]],
        derivative=c(1L,0L))
      h <- t(certificate$derivative.weights)
      attr(h,"operator.error") <- t(certificate$operator.error)
      return(h)
    }
    kw <- if(isTRUE(ingredients[[3L]])) NULL else
      .np_plot_lsq_gate_weights(bws,geometry,index,training=is.null(ex),package=package,max.bytes=max.bytes)
    .Call('C_np_plot_projection_apply',ingredients[[1L]],
      ingredients[[2L]][index,,drop=FALSE],kw,PACKAGE=package)
  }
}

# New-request GENERAL/SCALAR derivative adapter. ALL-LARGE positive-degree
# derivatives are a separate actual-owner construction, not qualified here.
.np_plot_lsq_derivative_factor <- function(bws, x, q, ex, coordinate,
    derivative.order = 1L, estimate, se, block.bytes = 4*1024^2,
    max.bytes = 128*1024^2) {
  n <- nrow(x); m <- nrow(ex)
  if (!is.data.frame(x) || !is.data.frame(ex) || n < 2L || m < 1L ||
      !identical(names(x), names(ex)) || length(q) != n ||
      !is.numeric(q) || any(!is.finite(q)) || anyNA(x) || anyNA(ex))
    stop('invalid LSQ derivative-factor data', call.=FALSE)
  continuous <- which(bws[['icon', exact=TRUE]])
  if(length(coordinate)!=1L || !is.finite(coordinate) ||
     !coordinate %in% continuous || length(derivative.order)!=1L ||
     !is.finite(derivative.order) || derivative.order<1L ||
     derivative.order!=as.integer(derivative.order))
    stop('invalid LSQ derivative request', call.=FALSE)
  if(length(estimate)!=m || length(se)!=m ||
     any(!is.finite(c(estimate,se))) || any(se<0))
    .np_plot_joint_unavailable('LSQ derivative marginal uncertainty is unavailable')
  if(length(block.bytes)!=1L || length(max.bytes)!=1L ||
     !is.finite(block.bytes) || !is.finite(max.bytes) ||
     block.bytes<=0 || max.bytes<=0)
    stop('invalid LSQ derivative-factor budget', call.=FALSE)
  spec <- npValidatedConditionalRegSpec(bws,
    where='LSQ actual derivative factor', ncon.field='ncon')
  scalar <- npIsCanonicalLp0Spec(spec, length(continuous))
  width <- if(scalar) 1 else npLpBasisNcol(spec$basis.engine,spec$degree.engine)
  fixed <- 64*1024^2+64*(as.double(n)+m)*(width+ncol(x)+1)+272*width^2+
    160*as.double(n)*m+64*as.double(n)
  if(identical(bws[['type', exact=TRUE]],'adaptive_nn'))
    fixed <- fixed+8*as.double(n)*min(n,64L)
  available <- min(block.bytes,max.bytes-fixed)
  rows <- min(n,m,floor(available/(160*as.double(n))))
  if(!is.finite(rows)||rows<1L)
    .np_plot_joint_unavailable('LSQ actual derivative factor exceeds the memory budget')
  # Native owner certification uses the complete evaluation geometry. Splitting
  # it merely to inspect ownership could change large-bandwidth eligibility.
  projection <- NULL
  if(!scalar) {
    projection <- .np_plot_lsq_gradient_projection(bws,x,q,ex,
      gradient.order=derivative.order)$fit.projection
    if(!is.list(projection)||length(projection)!=3L||
       length(projection[[3L]])!=1L||is.na(projection[[3L]]))
      stop('missing LSQ accepted owner certificate',call.=FALSE)
  }
  residual <- .np_plot_lsq_certified_residuals(bws,x,q,rows)
  z <- residual$residual
  derivative <- integer(length(continuous))
  derivative[match(coordinate,continuous)] <- as.integer(derivative.order)
  # Use the already implemented native derivative entry for higher orders too.
  # Public npreghat's generic higher-order dispatch is deliberately untouched.
  # Retain whole evaluation geometry; its full map is explicitly budgeted.
  D <- if(!is.null(projection) && isTRUE(projection[[3L]])) {
    # Public positive-degree fixed fits force tree OFF; ALL-LARGE is fixed-only.
    # This invariant certifies native donor identity without reimplementing it.
    if(!identical(bws[["type",exact=TRUE]],"fixed") ||
       !identical(.npreg_fit_tree_code(bws,bws$ncon,bws$nuno+bws$nord),DO_TREE_NO))
      stop("LSQ ALL-LARGE donor-order certificate is unavailable",call.=FALSE)
    .Call("C_np_plot_alllarge_derivative",projection[[1L]],
      as.matrix(x[,continuous,drop=FALSE]),as.matrix(ex[,continuous,drop=FALSE]),
      as.integer(spec$degree.engine),isTRUE(spec$bernstein.basis.engine),
      as.integer(npLpBasisCode(spec$basis.engine)),
      as.integer(match(coordinate,continuous)),as.integer(derivative.order),
      seq_len(n)-1L,projection[[2L]],PACKAGE="npRmpi")$weights
  } else if(!is.null(projection) && identical(bws[["ckertype",exact=TRUE]],"beta")) {
    certificate <- .np_plot_lsq_beta_moment_weights(bws,x,q,projection[[1L]],ex,
      max.bytes=max.bytes-fixed+64*as.double(n)*m,
      accepted.projection=projection[[2L]],
      derivative=as.integer(c(match(coordinate,continuous),derivative.order)))
    h <- t(certificate$derivative.weights)
    attr(h,"operator.error") <- t(certificate$operator.error)
    h
  } else .np_plot_lsq_gate_derivative(bws,x,ex,s=derivative,
    max.bytes=max.bytes-fixed+32*as.double(n)*m)
  if(!is.matrix(D)||!identical(dim(D),c(m,n))||any(!is.finite(D)))
    stop('invalid LSQ native derivative map',call.=FALSE)
  reconstructed <- drop(D%*%q)
  if(any(.np_plot_plreg_check_reconstruction(estimate,reconstructed,se)))
    .np_plot_joint_unavailable('LSQ gradient influence does not reconstruct the accepted gradient')
  A <- sweep(D,2L,z,'*')
  .np_plot_lsq_check_factor_precision(D,residual$interval,A)
  attr(A,"operator.error") <- NULL
  list(A=A,residual=z,estimate=reconstructed,block.rows=rows,
    forecast.bytes=fixed+160*n*rows,
    target='frozen LSQ transformed response; accepted GENERAL/SCALAR derivative map')
}

.np_plot_lsq_joint <- function(...) .npRmpi_with_local_regression(.np_plot_lsq_joint_local(...))

.np_plot_qreg_joint <- function(...) .npRmpi_with_local_regression(.np_plot_qreg_joint_local(...))

.np_plot_conmode_joint_columns <- function(...) .npRmpi_with_local_regression(.np_plot_conmode_joint_columns_local(...))

# New-inference only. The public LC categorical quantile target is
# -Delta F(q0,x)/f(q0,x), not the LP own-root quantile contrast.
.np_plot_qreg_lc_same_root_delta <- function(bws,x,y,ex,quantile,gradient,
                                    max.bytes=128*1024^2) {
  if (bws[['type',exact=TRUE]] %in% c('generalized_nn','adaptive_nn'))
    return(.np_plot_qreg_lc_same_root_nn_delta(bws,x,y,ex,quantile,gradient,max.bytes))
  spec <- npConditionalRegEngineSpec(bws,where='LC categorical quantile delta')
  if (!(identical(spec$reg.engine,'lc') || all(spec$degree.engine==0L)) ||
      length(gradient)!=1L || !is.finite(gradient) ||
      gradient!=as.integer(gradient) || gradient<1L || gradient>ncol(x) ||
      !is.factor(x[[gradient]]) || ncol(y)!=1L || !is.numeric(y[[1L]]))
    stop('invalid LC categorical quantile target',call.=FALSE)
  n<-nrow(x);m<-nrow(ex)
  if(n<2L || m<1L || nrow(y)!=n || length(quantile)!=m ||
     any(!is.finite(quantile)) || anyNA(x) || anyNA(y) || anyNA(ex))
    stop('invalid LC categorical quantile data',call.=FALSE)
  if(!identical(bws[['type',exact=TRUE]],'fixed') ||
     !identical(as.integer(bws[['cykerorder',exact=TRUE]]),2L) ||
     !bws[['cykertype',exact=TRUE]] %in% c('gaussian','epanechnikov','beta','uniform'))
    .np_plot_joint_unavailable('LC categorical full delta needs a qualified response derivative for this bandwidth/kernel geometry')
  forecast<-128*as.double(n)*m+64*(as.double(n)+m)*(ncol(x)+2)
  if(length(max.bytes)!=1L || !is.finite(max.bytes) || max.bytes<=0 ||
     !is.finite(forecast) || forecast>max.bytes)
    .np_plot_joint_unavailable('LC categorical full delta exceeds its memory budget')
  ey<-setNames(data.frame(quantile),names(y))
  F0<-.np_plot_density_conditional(bws,x,y,ex,ey,cdf=TRUE)
  D<-.np_plot_density_conditional(bws,x,y,ex,ey,cdf=TRUE,gradient=gradient)
  f0<-.np_plot_density_conditional(bws,x,y,ex,ey,cdf=FALSE)
  xb<-.npcdhat_make_xkbw(bws,x);yb<-.npcdhat_make_ybw(bws,y)
  K<-.np_plot_density_matrix(xb,x,ex,rep('normal',ncol(x)))
  kind<-yb[['ckertype',exact=TRUE]];h<-yb[['bw',exact=TRUE]][1L]
  lb<-yb[['ckerlb',exact=TRUE]][1L];ub<-yb[['ckerub',exact=TRUE]][1L]
  if(any(quantile<=lb | quantile>=ub))
    .np_plot_joint_unavailable('LC categorical quantile root is on a support boundary')
  # A compact density is nondifferentiable at an exact donor support face.
  # Interval-normalizer faces must be checked as well for bounded kernels.
  if(kind %in% c('epanechnikov','uniform')) {
    face<-if(kind=='uniform')1 else sqrt(5)
    z<-outer(quantile,y[[1L]],'-')/h
    norm.faces<-c(if(is.finite(lb))(quantile-lb)/h,
                  if(is.finite(ub))(quantile-ub)/h)
    if(any(abs(z)==face) || any(abs(norm.faces)==face))
      .np_plot_joint_unavailable('LC categorical quantile density is nondifferentiable at a compact support face')
  }
  # The CDF is donor-centered; the PDF is evaluation-centered. Swap the
  # canonical normal-kernel roles while retaining the fixed physical h/support.
  # This is essential for bounded and beta kernels; no symmetry is presumed.
  FY<-t(.np_plot_density_matrix(yb,ey,y,'normal'))
  frames<-npCategoricalFirstDifferenceFrames(ex,index=gradient,where='LC full delta')
  endpoint<-function(xx) {
    k<-.np_plot_density_matrix(xb,x,xx,rep('normal',ncol(x)))
    k<-.np_plot_density_factor_weights(k,xb,x,xx,rep('normal',ncol(x)),ordered=FALSE)
    rowSums(k*FY)/rowSums(k)
  }
  Dy<-endpoint(frames$upper)-endpoint(frames$lower)
  Fy<-rowSums(K*FY)/rowSums(K)
  dY<-.np_plot_density_derivative_matrix(yb,y,ey,1L)
  fy<-rowSums(K*dY)/rowSums(K)
  f<-f0$estimate
  # Translation-kernel OP_DERIVATIVE differentiates the numerator. The
  # evaluation-centered bounded PDF also has a moving truncation normalizer.
  # Canonical normalized boundary kernel values give exactly M'(q)/M(q).
  if(kind!='beta' && (is.finite(lb) || is.finite(ub))) {
    edge<-function(value) {
      if(!is.finite(value))return(rep(0,m))
      yy<-setNames(data.frame(value),names(y))
      drop(.np_plot_density_matrix(yb,yy,ey,'normal'))
    }
    fy<-fy-f*(edge(lb)-edge(ub))
  }
  if(any(!is.finite(c(F0$estimate,D$estimate,f,Fy,Dy,fy))) || any(f<=0 | Fy<=0))
    .np_plot_joint_unavailable('LC categorical full delta requires regular positive root and density slopes')
  Aq<--F0$A/Fy
  ratio<-D$estimate/f
  root.coefficient<- -Dy/f+ratio*(fy/f)
  A<- -D$A/f+ratio*(f0$A/f)+root.coefficient*Aq
  if(any(!is.finite(A)))
    .np_plot_joint_unavailable('LC categorical full delta is not representable')
  se<-vapply(seq_len(m),function(i){
    scale<-max(abs(A[i,]));if(scale==0)return(0)
    value<-scale*sqrt(sum((A[i,]/scale)^2))
    if(!is.finite(value)||value==0)
      .np_plot_joint_unavailable('LC categorical full-delta norm is not representable')
    value
  },numeric(1))
  list(A=A,se=se,estimate= -ratio,root=F0$estimate,
    primitives=list(AF=F0$A,AD=D$A,Af=f0$A,D=D$estimate,f=f,
      Fy=Fy,Dy=Dy,fy=fy,Aq=Aq),forecast.bytes=forecast,
    target='full same-central-root LC categorical quantile plug-in delta; fixed smoothing',
    variance='new simultaneous-request delta scale; incumbent marginal intervals unchanged')
}

# New-request marshaler only; no protected function is rebound or rewritten.
.np_plot_lsq_gradient_projection <- function(bws,x,q,ex,gradient.order=1L,
                                                package='npRmpi',return.args=FALSE) {
  x<-toFrame(x);ex<-toFrame(ex)
  if(!(x %~% ex) || nrow(x)<2L || nrow(ex)<1L || length(q)!=nrow(x) ||
     !is.numeric(q) || any(!is.finite(q)) || anyNA(x) || anyNA(ex) ||
     length(bws[['bw',exact=TRUE]])!=ncol(x))
    stop('invalid plot gradient projection data',call.=FALSE)
  spec<-npValidatedConditionalRegSpec(bws,where='plot gradient projection',ncon.field='ncon')
  if(bws$ncon<1L || npIsCanonicalLp0Spec(spec,bws$ncon))
    stop('plot gradient projection requires positive-degree LP',call.=FALSE)
  order<-npValidateGlpGradientOrder(regtype='lp',gradient.order=gradient.order,ncon=bws$ncon)
  if(!any(order<=spec$degree.engine))
    stop('plot gradient projection has no available derivative',call.=FALSE)
  reg<-npRegtypeToC(regtype=spec$regtype.engine,degree=spec$degree.engine,
    ncon=bws$ncon,context='plot gradient projection')
  beta<-identical(bws[['ckertype',exact=TRUE]],'beta')
  npValidateBetaKernelSpecification(ckertype=bws$ckertype,ckerorder=bws$ckerorder,
    bwtype=bws$type,ckerbound=bws$ckerbound,ckerlb=bws$ckerlb,ckerub=bws$ckerub,
    dati=bws$xdati,bw=bws$bw,bandwidth.compute=FALSE,where='plot gradient projection',
    regtype=bws$regtype,allow.categorical=TRUE,allow.general.lp=TRUE)
  x<-adjustLevels(x,bws$xdati);ex<-adjustLevels(ex,bws$xdati,allowNewCells=TRUE)
  npKernelBoundsCheckEval(ex,bws$icon,bws$ckerlb,bws$ckerub,argprefix='cker')
  xx<-toMatrix(x);ee<-toMatrix(ex)
  npCheckRegressionDesignCondition(reg.code=reg$code,xcon=xx[,bws$icon,drop=FALSE],
    basis=spec$basis.engine,degree=spec$degree.engine,
    bernstein.basis=spec$bernstein.basis.engine,where='plot gradient projection')
  opt<-list(num_obs_train=nrow(x),num_obs_eval=nrow(ex),num_uno=bws$nuno,
    num_ord=bws$nord,num_con=bws$ncon,int_LARGE_SF=SF_ARB,
    BANDWIDTH_reg_extern=switch(bws$type,fixed=BW_FIXED,generalized_nn=BW_GEN_NN,adaptive_nn=BW_ADAP_NN),
    int_MINIMIZE_IO=if(isTRUE(getOption('np.messages')))IO_MIN_FALSE else IO_MIN_TRUE,
    kerneval=switch(bws$ckertype,gaussian=CKER_GAUSS+bws$ckerorder/2-1,
      epanechnikov=CKER_EPAN+bws$ckerorder/2-1,uniform=CKER_UNI,beta=CKER_COORDINATE),
    ukerneval=switch(bws$ukertype,aitchisonaitken=UKER_AIT,liracine=UKER_LR),
    okerneval=switch(bws$okertype,wangvanryzin=OKER_WANG,liracine=OKER_LR,racineliyan=OKER_RLY),
    ey_is_ty=TRUE,do_grad=TRUE,regtype=reg$code,no.ex=FALSE,
    mcv.numRow=attr(bws$xmcv,'num.row'),
    int_do_tree=if(beta)DO_TREE_NO else .npreg_fit_tree_code(bws,bws$ncon,bws$nuno+bws$nord),
    categorical.compress=if(beta && bws$nuno+bws$nord>0L)
      npStrictLogicalOption('np.categorical.compress',TRUE) else FALSE)
  opt<-c(opt,npContinuousKernelDescriptorOptions(bws))
  bw<-.np_physical_bandwidth(bws)
  bounds<-npKernelBoundsMarshal(bws$ckerlb[bws$icon],bws$ckerub[bws$icon])
  args<-list(as.double(xx[,bws$iuno,drop=FALSE]),as.double(xx[,bws$iord,drop=FALSE]),
    as.double(xx[,bws$icon,drop=FALSE]),as.double(q),
    as.double(ee[,bws$iuno,drop=FALSE]),as.double(ee[,bws$iord,drop=FALSE]),
    as.double(ee[,bws$icon,drop=FALSE]),double(),
    as.double(c(bw[bws$icon],bw[bws$iuno],bw[bws$iord])),as.double(bws$xmcv),
    as.double(attr(bws$xmcv,'pad.num')),as.double(bws$nconfac),as.double(bws$ncatfac),
    as.double(bws$sdev),as.integer(opt),as.integer(reg$degree),as.integer(order),
    as.integer(isTRUE(spec$bernstein.basis.engine)),as.integer(npLpBasisCode(spec$basis.engine)),
    as.integer(nrow(ex)),as.integer(ncol(x)),
    as.integer(c(.np_regression_output_request(se=FALSE,gradients=TRUE),0L,
      npLpBasisNcol(spec$basis.engine,reg$degree))),as.double(bounds$lb),as.double(bounds$ub))
  if(isTRUE(return.args))return(args)
  out<-.np_with_nn_radius_context(
    .Call('C_np_plot_regression_gradient_projection',args,PACKAGE=package),
    continuous.names=bws$xnames[bws$icon],bws=bws,ntrain=nrow(x))
  reverse<-order(c(which(bws$icon),which(bws$iuno),which(bws$iord)))
  out$gradient<-matrix(out$gradient,nrow(ex),ncol(x))[,reverse,drop=FALSE]
  out
}

# New-inference integration only; caller must already own package local-MPI scope.
.np_plot_qreg_nn_response_jet <- function(bws,y,query,max.bytes=128*1024^2,package='npRmpi') {
  if(!is.data.frame(y)||ncol(y)!=1L||!is.numeric(y[[1L]])||
     !is.numeric(query)||any(!is.finite(query)))stop('invalid NN response frames')
  if(bws[['cykerorder',exact=TRUE]]!=2L)stop('public quantile response kernel order must be two')
  mode<-switch(bws[['type',exact=TRUE]],generalized_nn=1L,adaptive_nn=2L,
    stop('NN response jet requires an NN bandwidth'))
  kind<-switch(bws[['cykertype',exact=TRUE]],gaussian=0L,epanechnikov=4L,
    uniform=8L,beta=100L,stop('unknown NN quantile response kernel'))
  yb<-.npcdhat_make_ybw(bws,y)
  tryCatch(.Call('C_np_plot_nn_response_jet',as.double(y[[1L]]),as.double(query),
    as.double(yb[['bw',exact=TRUE]][1L]),mode,kind,
    as.double(c(yb[['ckerlb',exact=TRUE]][1L],yb[['ckerub',exact=TRUE]][1L])),
    as.double(max.bytes),PACKAGE=package),error=function(e) {
      message<-conditionMessage(e)
      for(prefix in c("np_plot_nn_unavailable: ","np_plot_beta_numerical_unavailable: "))
        if(startsWith(message,prefix))
          .np_plot_joint_unavailable(substring(message,nchar(prefix)+1L))
      stop(e)
    })
}

.np_plot_qreg_lc_same_root_nn_delta <- function(bws,x,y,ex,quantile,gradient,
                                      max.bytes=128*1024^2,package='npRmpi') {
  spec<-npConditionalRegEngineSpec(bws,where='LC NN categorical quantile delta')
  if(!(identical(spec$reg.engine,'lc')||all(spec$degree.engine==0L))||
     length(gradient)!=1L||!is.finite(gradient)||gradient!=as.integer(gradient)||
     gradient<1L||gradient>ncol(x)||!is.factor(x[[gradient]])||
     ncol(y)!=1L||!is.numeric(y[[1L]]))stop('invalid LC NN categorical target')
  # These positive bounds are numerical agreement/quadrature estimates, not
  # rigorous special-function enclosures. Do not lose nonzero discrepancies
  # through underflow before a reciprocal or derivative amplifies them.
  emul<-function(a,b) {
    z<-a*b
    if(any(!is.finite(z))||any(a>0&b>0&z==0))
      .np_plot_joint_unavailable("LC NN discrepancy multiplication is not representable")
    z
  }
  ediv<-function(a,b) {
    z<-a/b
    if(any(!is.finite(z))||any(a>0&z==0))
      .np_plot_joint_unavailable("LC NN discrepancy division is not representable")
    z
  }
  eadd<-function(a,b) {
    z<-a+b
    if(any(!is.finite(z))) .np_plot_joint_unavailable("LC NN discrepancy sum is not representable")
    z
  }
  esum<-function(a) {
    z<-rowSums(a);ops<-4*as.double(ncol(a))+4
    gamma<-ops*.Machine$double.eps/(1-ops*.Machine$double.eps)
    emul(z,1+gamma)
  }
  n<-nrow(x);m<-nrow(ex)
  forecast<-320*as.double(n)*m+128*(as.double(n)+m)*(ncol(x)+2)+12288
  if(n<2L||m<1L||nrow(y)!=n||length(quantile)!=m||anyNA(x)||anyNA(y)||anyNA(ex)||
     length(max.bytes)!=1L||!is.finite(max.bytes)||max.bytes<=0||
     !is.finite(forecast)||forecast>max.bytes)stop('invalid or oversized LC NN delta')
  ey<-setNames(data.frame(quantile),names(y))
  F0<-.np_plot_density_conditional(bws,x,y,ex,ey,cdf=TRUE)
  D<-.np_plot_density_conditional(bws,x,y,ex,ey,cdf=TRUE,gradient=gradient)
  f0<-.np_plot_density_conditional(bws,x,y,ex,ey,cdf=FALSE)
  jet<-.np_plot_qreg_nn_response_jet(bws,y,quantile,max.bytes,package)
  xb<-.npcdhat_make_xkbw(bws,x)
  K<-.np_plot_density_matrix(xb,x,ex,rep('normal',ncol(x)));den<-rowSums(K)
  Fy<-rowSums(K*jet$cdf.q)/den
  fy<-rowSums(K*jet$pdf.q)/den
  Fy.error<-ediv(esum(emul(abs(K),jet$cdf.q.error)),abs(den))
  f.error<-ediv(esum(emul(abs(K),jet$pdf.absolute.discrepancy)),abs(den))
  fy.error<-ediv(esum(emul(abs(K),jet$pdf.q.absolute.discrepancy)),abs(den))
  frames<-npCategoricalFirstDifferenceFrames(ex,index=gradient,where='LC NN full delta')
  endpoint<-function(xx) {
    k<-.np_plot_density_matrix(xb,x,xx,rep('normal',ncol(x)))
    k<-.np_plot_density_factor_weights(k,xb,x,xx,rep('normal',ncol(x)),ordered=FALSE)
    d<-rowSums(k)
    list(value=rowSums(k*jet$cdf.q)/d,error=ediv(esum(emul(abs(k),jet$cdf.q.error)),abs(d)))
  }
  u<-endpoint(frames$upper);l<-endpoint(frames$lower)
  Dy<-u$value-l$value;Dy.error<-eadd(u$error,l$error);f<-f0$estimate
  same.endpoint<-vapply(seq_len(m),function(j)identical(frames$upper[j,,drop=FALSE],
    frames$lower[j,,drop=FALSE]),logical(1))
  Dy.error[same.endpoint]<-0
  if(any(!is.finite(c(f,Fy,fy,Dy,Fy.error,Dy.error,f.error,fy.error)))||
     any(f<=f.error|Fy<=Fy.error))
    .np_plot_joint_unavailable('LC NN delta requires positive identifiable root/density slopes')
  if(any(abs(rowSums(K*jet$cdf)/den-F0$estimate)>2e-11*(1+abs(F0$estimate)))||
     any(abs(rowSums(K*jet$pdf)/den-f)>2e-11*(1+abs(f))))
    stop('NN response jet does not reconstruct the admitted conditional point map')
  raw.Af<-K*(jet$pdf-f)/den
  reconstructed.Af<-.np_plot_density_center(raw.Af)
  if(any(abs(reconstructed.Af-f0$A)>2e-11*(1+abs(f0$A))))
    stop('NN response jet does not reproduce the density influence owner')
  Af.error<-emul(ediv(abs(K),abs(den)),eadd(jet$pdf.absolute.discrepancy,f.error))
  Af.error<-eadd(Af.error,ediv(esum(Af.error),as.double(n)))
  Aq<--F0$A/Fy;ratio<-D$estimate/f
  coefficient<--Dy/f+ratio*fy/f
  A<--D$A/f+ratio*(f0$A/f)+coefficient*Aq
  if(any(!is.finite(A))) .np_plot_joint_unavailable('LC NN joint factor is not representable')
  se<-.np_plot_four_factor_se(A)
  # Propagate same-shape PDF amplitude disagreement and CDF concentration
  # quadrature through every full-delta term, including the random root and
  # density denominator. No existing point or influence amplitude is replaced.
  relative.f<-ediv(f.error,abs(f)-f.error)
  relative.Fy<-ediv(Fy.error,abs(Fy)-Fy.error)
  square.relative<-eadd(emul(2,relative.f),emul(relative.f,relative.f))
  reciprocal.square<-eadd(1,square.relative)
  Aqe<-emul(abs(Aq),relative.Fy)
  part1<-emul(abs(D$A/f),relative.f)
  part2<-emul(abs(ratio),eadd(emul(ediv(Af.error,abs(f)),reciprocal.square),
    emul(abs(f0$A/f),square.relative)))
  ce<-eadd(emul(ediv(Dy.error,abs(f)),eadd(1,relative.f)),
    emul(abs(Dy/f),relative.f))
  ce<-eadd(ce,emul(abs(ratio),eadd(emul(ediv(fy.error,abs(f)),reciprocal.square),
    emul(abs(fy/f),square.relative))))
  Ae<-eadd(eadd(part1,part2),eadd(emul(ce,eadd(abs(Aq),Aqe)),emul(abs(coefficient),Aqe)))
  factor.error<-.np_plot_four_factor_se(Ae)
  qualified.error<-ifelse(se>0,factor.error/se<=1e-8,factor.error==0)
  if(any(!qualified.error))
    .np_plot_joint_unavailable('LC NN numerical agreement or quadrature exceeds the joint-factor error budget')
  list(A=A,se=se,estimate=-ratio,root=F0$estimate,
    primitives=list(AF=F0$A,AD=D$A,Af=f0$A,D=D$estimate,f=f,Fy=Fy,Dy=Dy,fy=fy,Aq=Aq),
    jet=jet,numerical.factor.error=factor.error,numerical.factor.error.matrix=Ae,
    forecast.bytes=forecast,
    target='same-central-root LC categorical quantile full local delta; realized NN geometry, query radius motion retained',
    variance='conditional empirical-mass approximation; no unconditional selected-radius law asserted')
}

# Exact accepted-beta-moment scale replay.
# B must be the accepted fit's original-order training basis; q its exact response.
.np_plot_lsq_beta_moment_weights <- function(bws,x,q,B,ex=NULL,index=NULL,
                                  max.bytes=128*1024^2,package='npRmpi',accepted.projection=NULL,
                                  derivative=NULL) {
  training<-is.null(ex);x<-toFrame(x);if(training)ex<-x else ex<-toFrame(ex)
  if(!identical(bws[['ckertype',exact=TRUE]],'beta')||!is.matrix(B)||
     !is.double(B)||nrow(B)!=nrow(x)||ncol(B)<1L||!is.numeric(q)||length(q)!=nrow(x)||
     any(!is.finite(q))||!(x %~% ex)||anyNA(x)||anyNA(ex))stop('invalid beta moment replay inputs')
  if(is.null(index))index<-seq_len(nrow(ex))
  if(!is.integer(index)||!length(index)||anyNA(index)||any(index<1L|index>nrow(ex)))
    stop('invalid beta moment replay row map')
  n<-nrow(x);m<-length(index);w<-ncol(B);p<-ncol(x)
  forecast<-64*as.double(n)*m+128*(as.double(n)+m)*(p+w+3)+256*as.double(w)*(w+2)+16*as.double(w)*w+96*w+32
  if(length(max.bytes)!=1L||!is.finite(max.bytes)||max.bytes<=0||
     !is.finite(forecast)||forecast>max.bytes||as.double(n)*m>.Machine$integer.max)
    .np_plot_joint_unavailable('beta moment replay exceeds its resource budget')
  x<-adjustLevels(x,bws$xdati);ex<-adjustLevels(ex,bws$xdati,allowNewCells=TRUE)
  split<-function(dat){z<-toMatrix(dat);lapply(list(uno=bws$iuno,ord=bws$iord,con=bws$icon),
    function(i){a<-z[,i,drop=FALSE];storage.mode(a)<-'double';a})}
  bw<-.np_physical_bandwidth(bws)
  code<-as.integer(c(switch(bws$type,fixed=BW_FIXED,generalized_nn=BW_GEN_NN,adaptive_nn=BW_ADAP_NN),
    bws$ckerorder,switch(bws$ukertype,aitchisonaitken=UKER_AIT,liracine=UKER_LR),
    switch(bws$okertype,wangvanryzin=OKER_WANG,liracine=OKER_LR,nliracine=OKER_NLR,racineliyan=OKER_RLY),
    bws$nuno+bws$nord>0L&&npStrictLogicalOption('np.categorical.compress',TRUE)))
  args<-list(train=split(x),eval=split(ex[index,,drop=FALSE]),basis=B,response=as.double(q),
    bandwidth=as.double(c(bw[bws$icon],bw[bws$iuno],bw[bws$iord])),kernel=code,
    lower=as.double(bws$ckerlb[bws$icon]),upper=as.double(bws$ckerub[bws$icon]),
    support=.np_native_categorical_support(bws),eval.train.index=if(training)index-1L else integer(),
    max.bytes=as.double(max.bytes))
  if(!is.null(derivative)){
    spec<-npValidatedConditionalRegSpec(bws,where='beta derivative replay',ncon.field='ncon')
    if(!is.integer(derivative)||length(derivative)!=2L||anyNA(derivative)||
       derivative[1]<1L||derivative[1]>bws$ncon||derivative[2]<0L||
       !is.matrix(accepted.projection)||!is.double(accepted.projection)||
       !identical(dim(accepted.projection),c(nrow(ex),w)))
      stop('invalid beta derivative replay certificate')
    args<-c(args,list(accepted.projection=accepted.projection[index,,drop=FALSE],
      degree=as.integer(spec$degree.engine),basis.code=as.integer(npLpBasisCode(spec$basis.engine)),
      bernstein=as.integer(isTRUE(spec$bernstein.basis.engine)),
      coordinate=derivative[1]-1L,order=derivative[2]))
  }
  result<-.np_with_nn_radius_context(.Call('C_np_plot_beta_moment_weights',args,PACKAGE=package),
    continuous.names=bws$xnames[bws$icon],bws=bws,ntrain=n)
  if(!is.null(derivative)&&any(!result[['certified',exact=TRUE]]))
    .np_plot_joint_unavailable('accepted beta operator cannot certify requested floating-point precision')
  result
}


.np_plot_ann_rank_geometry <- function(bws,data,package='npRmpi') {
  x<-as.matrix(data[,bws[['icon',exact=TRUE]],drop=FALSE]);storage.mode(x)<-'double'
  raw<-as.double(unlist(bws[['bandwidth',exact=TRUE]],use.names=FALSE))[bws[['icon',exact=TRUE]]]
  rank<-.Call('C_np_plot_ann_rank_decode',raw,as.integer(nrow(x)),PACKAGE=package)
  if(any(rank$extended)||any(rank$distance.scale!=1))
    stop('np_plot_ann_unavailable: extended maximum-radius Gaussian uncertainty requires a separate endpoint or tail law')
  orders<-lapply(seq_len(ncol(x)),function(j)order(x[,j]))
  geometry<-lapply(seq_len(ncol(x)),function(j).Call('C_np_ann_geometry',x[orders[[j]],j],rank$lookup[j],PACKAGE=package))
  list(x=x,order=orders,geometry=geometry,lookup=rank$lookup,distance.scale=rank$distance.scale,
    requested=raw,valid=all(vapply(geometry,function(g)all(g[['status']]==0L),logical(1))))
}

.np_plot_ann_bounded_density_factor <- function(bws,data,evaluation,owner,
  package='npRmpi',max.bytes=128*1024^2,tile.rows=NULL) {
  stopifnot(is.environment(owner));f<-function(s)get(s,owner,inherits=TRUE)
  if(!identical(bws[['type',exact=TRUE]],'adaptive_nn')||
     !is.data.frame(data)||!is.data.frame(evaluation)||anyNA(data)||anyNA(evaluation)||
     !identical(names(data),bws[['xnames',exact=TRUE]])||
     !identical(names(evaluation),names(data)))stop('invalid ANN plot data contract')
  positive.mul<-function(a,b) {
    z<-a*b
    if(any(!is.finite(z))||any(a>0&b>0&z==0))stop('np_plot_beta_numerical_unavailable: discrepancy multiplication not representable')
    z
  }
  positive.div<-function(a,b) {
    z<-a/b
    if(any(!is.finite(z))||any(a>0&z==0))stop('np_plot_beta_numerical_unavailable: discrepancy division not representable')
    z
  }
  n<-nrow(data);m<-nrow(evaluation);p<-ncol(data);con<-which(bws[['icon',exact=TRUE]]);pc<-length(con)
  if(n<2L||m<1L||pc<1L||length(max.bytes)!=1L||!is.finite(max.bytes)||max.bytes<=0)
    stop('invalid ANN plot factor dimensions or budget')
  tile<-min(m,64,max(1,floor(4*1024^2/(8*as.double(n)*(20+8*p)))))
  if(!is.null(tile.rows)) {
    if(length(tile.rows)!=1L||!is.finite(tile.rows)||tile.rows!=floor(tile.rows)||tile.rows<1||tile.rows>tile)
      stop('tile.rows may only reduce the bounded tile')
    tile<-tile.rows
  }
  forecast<-48*as.double(n)*m+64*1024^2+8*as.double(n)*tile*(20+8*p)+
    256*as.double(n)*(p+1)+64*as.double(m)*(p+pc+1)
  if(!is.finite(forecast)||forecast>max.bytes)stop('np_plot_ann_unavailable: ANN plot factor exceeds resident memory budget')
  # Validate retained/physical metadata, then invoke the canonical native
  # decoder before the unchanged integer-rank spacing pilot.
  raw<-as.double(unlist(bws[['bandwidth',exact=TRUE]],use.names=FALSE))
  physical<-as.double(f('.np_physical_bandwidth')(bws))
  if(length(raw)!=p||length(physical)!=p||any(!is.finite(raw))||
     !identical(raw[con],physical[con])||any(raw[con]<1))stop('invalid canonical ANN rank metadata')
  st<-.np_plot_ann_rank_geometry(bws,data,package)
  if(!isTRUE(st$valid))stop('np_plot_ann_unavailable: ANN plot radial spacing geometry is nonregular')
  kb<-f('.np_make_kbandwidth_unconditional')(bws,data)
  order<-as.integer(bws[['ckerorder',exact=TRUE]])
  ker<-switch(bws[['ckertype',exact=TRUE]],gaussian=match(order,c(2L,4L,6L,8L))-1L,
    epanechnikov=match(order,c(2L,4L,6L,8L))+3L,uniform=8L,beta=match(order,c(2L,4L,6L,8L))+99L,stop('unknown ANN kernel'))
  if(length(ker)!=1L||is.na(ker))stop('unknown ANN kernel order')
  x<-as.matrix(data[con]);storage.mode(x)<-'double'
  h<-matrix(NA_real_,n,pc)
  for(j in seq_len(pc))h[st$order[[j]],j]<-st$geometry[[j]]$radius
  # Fitted kernel bounds retain original-variable indexing, including factors.
  bounds<-rbind(as.double(kb[['ckerlb',exact=TRUE]][con]),as.double(kb[['ckerub',exact=TRUE]][con]))
  if(!identical(dim(bounds),c(2L,as.integer(pc))))stop('native continuous bounds order is malformed')
  mb<-lapply(seq_len(p),function(j)f('.npcopula_marginal_bw')(bws,data,j,target='density',kbandwidth=TRUE))
  train<-lapply(seq_len(p),function(j)f('.npcopula_marginal_data')(bws,data,j))
  evaluate<-lapply(seq_len(p),function(j)f('.npcopula_marginal_eval_data')(bws,evaluation,j))
  A<-error.A<-matrix(NA_real_,m,n);estimate<-denominator<-numeric(m);face.record<-list()
  for(start in seq.int(1L,m,by=tile)) {
    rows<-seq.int(start,min(m,start+tile-1L));nr<-length(rows)
    ev<-as.matrix(evaluation[rows,con,drop=FALSE]);storage.mode(ev)<-'double'
    jet<-.Call('C_np_plot_ann_kernel_terms',x,ev,h,as.integer(ker),bounds,as.double(max.bytes),PACKAGE=package)
    W<-WE<-vector('list',p)
    for(j in seq_len(pc)) {
      W[[con[j]]]<-matrix(jet$weight[,,j],n,nr)
      WE[[con[j]]]<-matrix(jet$weight.absolute.discrepancy[,,j],n,nr)
    }
    for(j in setdiff(seq_len(p),con)) {
      z<-f('.np_plot_density_matrix')(mb[[j]],train[[j]],evaluate[[j]][rows,,drop=FALSE],'normal')
      W[[j]]<-t(f('.np_plot_density_factor_weights')(z,mb[[j]],train[[j]],evaluate[[j]][rows,,drop=FALSE],'normal'));WE[[j]]<-matrix(0,n,nr)
    }
    product<-function(keep)if(length(keep))Reduce('*',W[keep])else matrix(1,n,nr)
    product.error<-function(keep) {
      value<-matrix(1,n,nr);error<-matrix(0,n,nr)
      for(j in keep) {error<-positive.mul(error,abs(W[[j]]))+positive.mul(abs(value),WE[[j]])+positive.mul(error,WE[[j]]);value<-value*W[[j]]}
      error
    }
    K<-product(seq_len(p));KE<-product.error(seq_len(p));D<-rep(1,nr);DE<-rep(0,nr)
    for(j in seq_len(p)) {
      dj<-f('.np_ksum_unconditional_eval_exact')(train[[j]],evaluate[[j]][rows,,drop=FALSE],mb[[j]],'normal')
      if(length(dj)!=nr)stop('ANN marginal row count mismatch')
      if(any(!is.finite(dj))||any(dj<=0))stop('np_plot_ann_unavailable: copula needs positive raw marginal densities')
      de<-positive.div(colSums(WE[[j]]),n)
      DE<-positive.mul(DE,abs(dj))+positive.mul(abs(D),de)+positive.mul(DE,de)
      D<-D*dj
      if(any(!is.finite(D))||any(D<=0))stop('np_plot_ann_unavailable: ANN copula denominator is not representable and positive')
    }
    radius<-lapply(seq_len(pc),function(j)matrix(jet$negative.h.derivative[,,j],n,nr)*product(setdiff(seq_len(p),con[j])))
    radius.error<-lapply(seq_len(pc),function(j) {
      keep<-setdiff(seq_len(p),con[j]);v<-matrix(jet$negative.h.derivative[,,j],n,nr);e<-matrix(jet$radius.absolute.discrepancy[,,j],n,nr)
      other<-product(keep);oe<-product.error(keep)
      positive.mul(e,abs(other))+positive.mul(abs(v),oe)+positive.mul(e,oe)
    })
    if(any(!is.finite(DE)))stop('np_plot_beta_numerical_unavailable: marginal amplitude discrepancy is not representable')
    if(any(DE<0))stop('negative ANN marginal discrepancy bound')
    if(any(DE>=D))stop('np_plot_beta_numerical_unavailable: marginal amplitude agreement failed')
    faces<-NULL
    if(ker==8L) {
      cut<-coefficient<-matrix(NA_real_,nr,2L*pc)
      for(j in seq_len(pc)) {
        z<-st$x[st$order[[j]],j];marks<-product(setdiff(seq_len(p),con[j]))
        rank<-as.double(st$lookup[j])
        for(q in seq_len(nr))for(side in 1:2) {
          column<-2L*(j-1L)+side;xq<-ev[q,j]
          probability<-findInterval(xq,z)/n+c(-1,1)[side]*rank/(n-1)
          if(probability<0||probability>1) {cut[q,column]<-xq;coefficient[q,column]<-0;next}
          if(probability==0||probability==1)stop('np_plot_ann_unavailable: ANN uniform face is at a probability endpoint')
          at<-1+(n-1)*probability;lo<-floor(at);s<-z[lo]+(at-lo)*(z[lo+1L]-z[lo])
          center<-(xq+s)/2;r<-abs(xq-s)/2
          fs<-f('.np_ann_face_window')(z,s);gt<-f('.np_ann_face_window')(z,center)
          if(is.null(fs)||is.null(gt)||!is.finite(r)||r<=0)stop('np_plot_ann_unavailable: ANN uniform face pilot is nonregular')
          ids<-st$order[[j]][seq.int(gt$lo,gt$hi)];mark<-marks[ids,q]
          mark[c(1L,length(mark))]<-mark[c(1L,length(mark))]/2
          marked<-(sum(mark)/n)/gt$width
          mass<-.Call('C_np_plot_ann_kernel_terms',matrix(center,1,1),matrix(numeric(),0,1),
            matrix(r,1,1),8L,matrix(bounds[,j],2,1),as.double(max.bytes),PACKAGE=package)$mass[1L]
          value<-c(-1,1)[side]*(.25*marked/fs$density)/r/mass
          if(!is.finite(value))stop('np_plot_ann_unavailable: ANN uniform face coefficient is not representable')
          cut[q,column]<-s;coefficient[q,column]<-value
        }
      }
      faces<-list(cut=cut,coefficient=coefficient)
      face.record[[as.character(start)]]<-faces
    }
    for(q in seq_len(nr)) {
      phi<-K[,q];ep<-KE[,q]
      for(j in seq_len(pc)) {
        g<-st$geometry[[j]];o<-st$order[[j]];w<-radius[[j]][o,q]/(g$radius*g$density.sum)
        delta<-numeric(n+1L);ag<-rowsum(c(w,-w),c(g$left,g$right)+1L,reorder=FALSE)
        delta[as.integer(rownames(ag))]<-ag[,1L];phi[o]<-phi[o]+cumsum(delta)[seq_len(n)]/n
        we<-positive.div(radius.error[[j]][o,q],g$radius*g$density.sum)
        de<-numeric(n+1L);ae<-rowsum(c(we,-we),c(g$left,g$right)+1L,reorder=FALSE)
        de[as.integer(rownames(ae))]<-ae[,1L]
        total.e<-cumsum(de)[seq_len(n)]
        ee<-positive.div(abs(total.e),n)*sign(total.e)
        # Roundoff allowance applies only to a nonnegative uncertainty sum.
        ops<-4*as.double(n)+4;gamma<-ops*.Machine$double.eps/(1-ops*.Machine$double.eps)
        allowance<-positive.div(positive.mul(gamma,sum(we)),n)
        if(any(ee < -allowance))stop('negative ANN amplitude discrepancy accumulation')
        ep[o]<-ep[o]+pmax(0,ee)+allowance
        if(!is.null(faces))for(side in 1:2) {
          column<-2L*(j-1L)+side
          phi<-phi+faces$coefficient[q,column]*((st$x[,j]<=ev[q,j])-(st$x[,j]<=faces$cut[q,column]))
        }
      }
      unscaled<-(phi-mean(phi))/n*sqrt(n/(n-1))
      ue<-positive.mul(positive.div(ep+positive.div(sum(ep),n),n),sqrt(n/(n-1)))
      A[rows[q],]<-unscaled/D[q]
      error.A[rows[q],]<-positive.div(ue,D[q]-DE[q])+positive.div(positive.mul(abs(unscaled),DE[q]),D[q]*(D[q]-DE[q]))
    }
    estimate[rows]<-colMeans(K)/D;denominator[rows]<-D
  }
  if(any(!is.finite(c(A,estimate,denominator))))stop('np_plot_ann_unavailable: ANN copula factor is not representable')
  scale<-f('.np_plot_four_factor_se')(A);bound<-f('.np_plot_four_factor_se')(error.A)
  ratio<-ifelse(scale>0,bound/scale,ifelse(bound==0,0,Inf))
  if(any(!is.finite(ratio))||any(ratio>1e-8))stop('np_plot_beta_numerical_unavailable: factor-relative amplitude disagreement')
  list(A=A,estimate=estimate,denominator=denominator,
    simultaneous.se=scale,exact.diagonal=TRUE,
    factor.absolute.discrepancy=error.A,factor.discrepancy.norm=bound,factor.discrepancy.ratio=ratio,
    numerical.contract='same-shape canonical PDF versus independent Rmath agreement; not rigorous special-function enclosure',
    target='ANN copula density; frozen inverse coordinates and positive marginal denominator; normalized numerator random-radius working variability',
    marginal.contract='new simultaneous scale only; inherited se/pointwise/Bonferroni unchanged',
    face=face.record,forecast.bytes=forecast,tile.rows=tile)
}

.np_plot_ann_bounded_distribution_factor <- function(bws,data,evaluation,owner,
  package='npRmpi',max.bytes=128*1024^2,tile.rows=NULL) {
  stopifnot(is.environment(owner));f<-function(s)get(s,owner,inherits=TRUE)
  if(!identical(bws[['type',exact=TRUE]],'adaptive_nn')||
     !is.data.frame(data)||!is.data.frame(evaluation)||anyNA(data)||anyNA(evaluation)||
     !identical(names(data),bws[['xnames',exact=TRUE]])||
     !identical(names(evaluation),names(data)))stop('invalid ANN plot data contract')
  positive.mul<-function(a,b) {
    z<-a*b
    if(any(!is.finite(z))||any(a>0&b>0&z==0))stop('np_plot_cdf_numerical_unavailable: quadrature uncertainty multiplication not representable')
    z
  }
  positive.div<-function(a,b) {
    z<-a/b
    if(any(!is.finite(z))||any(a>0&z==0))stop('np_plot_cdf_numerical_unavailable: quadrature uncertainty division not representable')
    z
  }
  n<-nrow(data);m<-nrow(evaluation);p<-ncol(data);con<-which(bws[['icon',exact=TRUE]]);pc<-length(con)
  if(n<2L||m<1L||pc<1L||length(max.bytes)!=1L||!is.finite(max.bytes)||max.bytes<=0)
    stop('invalid ANN plot factor dimensions or budget')
  tile<-min(m,64,max(1,floor(4*1024^2/(8*as.double(n)*(16+6*p)))))
  if(!is.null(tile.rows)) {
    if(length(tile.rows)!=1L||!is.finite(tile.rows)||tile.rows!=floor(tile.rows)||tile.rows<1||tile.rows>tile)
      stop('tile.rows may only reduce the bounded tile')
    tile<-tile.rows
  }
  forecast<-32*as.double(n)*m+64*1024^2+8*as.double(n)*tile*(16+6*p)+
    256*as.double(n)*(p+1)+64*as.double(m)*(p+pc+1)
  if(!is.finite(forecast)||forecast>max.bytes)stop('np_plot_ann_unavailable: ANN plot factor exceeds resident memory budget')
  # Validate retained/physical metadata, then invoke the canonical native
  # decoder before the unchanged integer-rank spacing pilot.
  raw<-as.double(unlist(bws[['bandwidth',exact=TRUE]],use.names=FALSE))
  physical<-as.double(f('.np_physical_bandwidth')(bws))
  if(length(raw)!=p||length(physical)!=p||any(!is.finite(raw))||
     !identical(raw[con],physical[con])||any(raw[con]<1))stop('invalid canonical ANN rank metadata')
  st<-.np_plot_ann_rank_geometry(bws,data,package)
  if(!isTRUE(st$valid))stop('np_plot_ann_unavailable: ANN plot radial spacing geometry is nonregular')
  kb<-f('.np_make_kbandwidth_unconditional')(bws,data)
  order<-as.integer(bws[['ckerorder',exact=TRUE]])
  ker<-switch(bws[['ckertype',exact=TRUE]],gaussian=match(order,c(2L,4L,6L,8L))-1L,
    epanechnikov=match(order,c(2L,4L,6L,8L))+3L,uniform=8L,beta=match(order,c(2L,4L,6L,8L))+99L,stop('unknown ANN kernel'))
  if(length(ker)!=1L||is.na(ker))stop('unknown ANN kernel order')
  x<-as.matrix(data[con]);storage.mode(x)<-'double'
  h<-matrix(NA_real_,n,pc)
  for(j in seq_len(pc))h[st$order[[j]],j]<-st$geometry[[j]]$radius
  # Fitted kernel bounds retain original-variable indexing, including factors.
  bounds<-rbind(as.double(kb[['ckerlb',exact=TRUE]][con]),as.double(kb[['ckerub',exact=TRUE]][con]))
  if(!identical(dim(bounds),c(2L,as.integer(pc))))stop('native continuous bounds order is malformed')
  mb<-lapply(seq_len(p),function(j)f('.npcopula_marginal_bw')(bws,data,j,target='distribution',kbandwidth=TRUE))
  train<-lapply(seq_len(p),function(j)f('.npcopula_marginal_data')(bws,data,j))
  evaluate<-lapply(seq_len(p),function(j)f('.npcopula_marginal_eval_data')(bws,evaluation,j))
  A<-error.A<-matrix(NA_real_,m,n);estimate<-numeric(m)
  for(start in seq.int(1L,m,by=tile)) {
    rows<-seq.int(start,min(m,start+tile-1L));nr<-length(rows)
    ev<-as.matrix(evaluation[rows,con,drop=FALSE]);storage.mode(ev)<-'double'
    jet<-.Call('C_np_plot_ann_cdf_terms',x,ev,h,as.integer(ker),bounds,as.double(max.bytes),PACKAGE=package)
    W<-vector('list',p)
    for(j in seq_len(pc))W[[con[j]]]<-matrix(jet$cdf[,,j],n,nr)
    for(j in setdiff(seq_len(p),con)) {
      z<-f('.np_plot_density_matrix')(mb[[j]],train[[j]],evaluate[[j]][rows,,drop=FALSE],'integral')
      W[[j]]<-t(f('.np_plot_density_factor_weights')(z,mb[[j]],train[[j]],evaluate[[j]][rows,,drop=FALSE],'integral'))
    }
    product<-function(keep)if(length(keep))Reduce('*',W[keep])else matrix(1,n,nr)
    K<-product(seq_len(p))
    radius<-lapply(seq_len(pc),function(j)matrix(jet$negative.h.derivative[,,j],n,nr)*product(setdiff(seq_len(p),con[j])))
    radius.error<-lapply(seq_len(pc),function(j)positive.mul(matrix(jet$radius.absolute.error[,,j],n,nr),abs(product(setdiff(seq_len(p),con[j])))))
    for(q in seq_len(nr)) {
      phi<-K[,q];ep<-numeric(n)
      for(j in seq_len(pc)) {
        g<-st$geometry[[j]];o<-st$order[[j]]
        divisor<-g$radius*g$density.sum
        w<-radius[[j]][o,q]/divisor;we<-positive.div(radius.error[[j]][o,q],divisor)
        prefix<-function(w) {
          delta<-numeric(n+1L);ag<-rowsum(c(w,-w),c(g$left,g$right)+1L,reorder=FALSE)
          delta[as.integer(rownames(ag))]<-ag[,1L];cumsum(delta)[seq_len(n)]/n
        }
        phi[o]<-phi[o]+prefix(w)
        # Every interval incidence is nonnegative. Prefix endpoint cancellation
        # computes the sum of absolute donor uncertainty, not a signed bound.
        delta.e<-numeric(n+1L);ae<-rowsum(c(we,-we),c(g$left,g$right)+1L,reorder=FALSE)
        delta.e[as.integer(rownames(ae))]<-ae[,1L]
        total.e<-cumsum(delta.e)[seq_len(n)]
        contribution<-positive.div(abs(total.e),n)*sign(total.e)
        ops<-4*as.double(n)+4;gamma<-ops*.Machine$double.eps/(1-ops*.Machine$double.eps)
        allowance<-positive.div(positive.mul(gamma,sum(we)),n)
        if(any(contribution< -allowance))stop('negative ANN quadrature bound accumulation')
        ep[o]<-ep[o]+pmax(0,contribution)+allowance
      }
      scale<-sqrt(n/(n-1))/n
      A[rows[q],]<-(phi-mean(phi))*scale
      error.A[rows[q],]<-positive.mul(ep+positive.div(sum(ep),n),scale)
    }
    estimate[rows]<-colMeans(K)
  }
  if(any(!is.finite(c(A,error.A,estimate))))stop('np_plot_ann_unavailable: ANN CDF factor/error is not representable')
  if(any(error.A<0))stop('negative ANN CDF factor error bound')
  scale<-f('.np_plot_four_factor_se')(A);bound<-f('.np_plot_four_factor_se')(error.A)
  ratio<-ifelse(scale>0,bound/scale,ifelse(bound==0,0,Inf))
  if(any(!is.finite(ratio))||any(ratio>1e-8))stop('np_plot_ann_unavailable: ANN CDF propagated quadrature error exceeds factor-relative bound')
  list(A=A,estimate=estimate,simultaneous.se=scale,exact.diagonal=TRUE,
    factor.absolute.error=error.A,factor.error.norm=bound,factor.error.ratio=ratio,
    target='ANN copula distribution; frozen inverse coordinates; same-sample random-radius numerator variability',
    marginal.contract='new simultaneous scale only; inherited se/pointwise/Bonferroni unchanged',
    uniform.contract='canonical open-K at isolated empirical CDF knots; no density face term; not finite-map q derivative',
    forecast.bytes=forecast,tile.rows=tile)
}
