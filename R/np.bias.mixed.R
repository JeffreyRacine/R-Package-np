# Frozen M1: convex mixtures of nested weighted polynomial projections.
# Thin bases and bounded observation rows only; no training-square pilot.
.np_bias_mixed_policy <- function(spec, bws, law) {
  mode <- spec$settings$adaptation_request
  requested <- spec$settings$requested
  defaults <- c(p0=4,p1=1.5,operator=4/3)
  override <- any(vapply(names(defaults),function(role) {
    key <- paste0(role,".scale"); value <- requested$bandwidth.controls[[key]]
    key %in% requested$supplied || is.null(value) ||
      !identical(as.double(value),as.double(defaults[[role]])) ||
      !is.null(names(value)) || length(value)!=1L
  },logical(1)))
  eligible <- length(spec$factor.cols)>0L && bws[["type"]]=="fixed" &&
    bws[["ckertype"]]=="gaussian" && bws[["ckerorder"]]==2L &&
    bws[["ckerbound"]]=="none" && isTRUE(spec$bernstein) && ncol(spec$Q)<=384L &&
    (is.null(law) || law %in% c("gaussian","mammen"))
  reason <- if(mode=="none") "disabled" else if(!eligible) "incumbent model or calibration law" else
    if(mode=="auto" && override) "explicit auxiliary geometry override" else "eligible"
  if(mode=="local-to-global" && !eligible)
    stop("mixed local-to-global adaptation requires fixed Gaussian GLP, stable auxiliary basis of at most 384 columns and Gaussian or Mammen calibration",call.=FALSE)
  active <- FALSE
  if(reason=="eligible") {
    b0 <- if(isTRUE(spec$pooled)) spec$h$p0 else spec$b
    active <- any(unname(bws[["bandwidth"]][["x"]][spec$numeric.cols])/spec$scale>b0)
    if(!active) reason <- "local geometry unchanged"
  }
  list(requested=mode,eligible=eligible,active=active,reason=reason,
    version="M1-nested-projection-1")
}

.np_bias_mixed_prepare <- function(s, bws, x, policy) {
  pooled <- isTRUE(s$pooled); q <- if(pooled)s$degree else s$q
  p <- npValidatedConditionalRegSpec(bws,where="mixed adaptation",ncon.field="ncon")[["degree.engine"]]
  geo <- .np_bias_adaptation_geometry(p,
    unname(bws[["bandwidth"]][["x"]][s$numeric.cols])/s$scale,
    if(pooled)s$h$p0 else s$b,if(pooled)s$h$p1 else s$p1$b,
    if(pooled)s$h$operator else s$h,q)
  groups <- if(pooled)list(seq_len(nrow(s$x))) else s$groups
  # The payload is response independent and can be serialized to MPI workers.
  blocks <- lapply(groups,function(ii) {
    Q <- s$Q[ii,,drop=FALSE]
    lambda <- c(geo$lambda,rep(0,ncol(Q)-length(geo$lambda)))
    list(indices=ii,Q=Q,x=s$x[ii,,drop=FALSE],xcon=as.matrix(s$xcon[ii,,drop=FALSE]),
      keys=s$keys[ii],pooled=pooled,
      spec=if(pooled)s[c("xcon","degree","main.degree","nc","fc","levels","bernstein")] else NULL,
      reference.x=s$xcon,nc=s$numeric.cols,fc=s$factor.cols,d=s$d,q=q,
      geo=geo,lambda=lambda,order=order(lambda,seq_along(lambda)),levels=sort(unique(lambda)))
  })
  s$mixed <- list(x=s$x,blocks=blocks,geometry=geo,iterations=s$settings$iterations,
    version=policy$version)
  owner <- .np_bias_mixed_p0_owner(s$mixed,progress=TRUE)
  z <- .np_bias_residual_apply(owner,matrix(s$y,ncol=1L))
  s$residual_raw <- drop(z$raw); s$residual <- drop(z$residual)
  s$residual_row_variance <- z$row_variance; s$mu <- s$y-s$residual_raw
  smoother <- .np_bias_student_variance_owner(list(x=s$xcon,groups=s$groups,h=s$variance_h))
  s$variance <- numeric(length(s$y))
  for(ii in .np_bias_profile_tiles(length(s$y)))
    s$variance[ii] <- drop(smoother(ii)%*%s$residual^2)
  if(any(!is.finite(s$variance)) || any(s$variance<=0))
    stop("nonpositive mixed adapted variance",call.=FALSE)
  s$main.request <- list(bws=bws,xdat=x,ydat=s$y)
  s$settings$policy <- paste0(s$settings$policy,"+M1-nested-projection-v1")
  s$settings$reference_bandwidth <- s$settings$bandwidth
  s$settings$reference_package_scale <- s$settings$effective_package_scale
  effective <- list(p0=geo$b0,p1=geo$b1,operator=geo$bK,variance=s$variance_h)
  for(role in names(effective)) {
    s$settings$bandwidth[[role]] <- setNames(effective[[role]]*s$scale,names(s$xcon))
    s$settings$effective_package_scale[[role]] <-
      s$settings$bandwidth[[role]]/(s$scale*nrow(x)^(-s$settings$exponent))
  }
  policy$geometry <- geo; policy$residual_policy <- "adapted P0 row-variance normalized"
  s$settings$adaptation <- policy
  s
}

.np_bias_mixed_basis <- function(b, at, axis=0L, order=1L) {
  nu <- integer(b$d); if(axis)nu[axis] <- order
  if(b$pooled) .np_bias_profile_basis(b$spec,at,if(axis)nu else NULL)
  else W.lp(b$reference.x,exdat=at[b$nc],degree=b$q,
    gradient.vec=if(axis)nu else NULL,basis="glp",bernstein.basis=TRUE)
}

.np_bias_mixed_kernel <- function(b, at, h, category=FALSE) {
  K <- .np_bias_adaptation_kernel(b$xcon,as.matrix(at[b$nc]),h)
  if(category && b$pooled)
    K <- K*outer(b$keys,.np_bias_profile_keys(at,b$fc),"==")
  K
}

.np_bias_mixed_pilot <- function(b, at, h, axis=0L, order=1L) {
  if(nrow(at)!=1L)stop("mixed pilot requires one anchor",call.=FALSE)
  g <- drop(.np_bias_mixed_kernel(b,at,h)); Q <- b$Q; P <- ncol(Q)
  value <- .np_bias_mixed_basis(b,at,axis,order)
  if(axis) {
    # Differentiate the complete projection, including its kernel weights.
    out <- numeric(nrow(Q)); lev <- b$levels; d <- 1/(1+lev); mix <- d-c(d[-1L],0)
    for(k in seq_along(lev)) {
      jj <- which(b$lambda<=lev[k])
      values <- do.call(rbind,lapply(0:order,function(r)
        .np_bias_mixed_basis(b,at,if(r)axis else 0L,r)[,jj,drop=FALSE]))
      ratios <- .np_bias_gaussian_derivatives((b$xcon[,axis]-at[[b$nc[axis]]])/h[axis],h[axis],order)
      row <- .np_bias_weighted_derivative(Q[,jj,drop=FALSE],g,values,ratios,
        role="M1 subspace derivative",rank.tolerance=1e-11)
      out <- out+mix[k]*row
    }
    return(list(H=matrix(out,1L,nrow(Q)),g=g))
  }
  jj <- b$order; Z <- Q[,jj,drop=FALSE]*sqrt(g); norms <- sqrt(colSums(Z^2))
  if(any(!is.finite(norms)) || any(norms<=0))
    stop("M1 nested subspace has an unsupported column",call.=FALSE)
  fit <- qr(sweep(Z,2L,norms,"/"),tol=1e-11,LAPACK=FALSE)
  if(fit$rank!=P || !identical(fit$pivot,seq_len(P)))
    stop("M1 nested subspace is numerically unidentified",call.=FALSE)
  U <- qr.Q(fit); R <- qr.R(fit); d <- 1/(1+b$lambda[jj])
  A <- backsolve(R,t(U)*d)*rep(sqrt(g),each=P)
  A <- sweep(A,1L,norms,"/"); coefficient <- A; coefficient[jj,] <- A
  list(A=coefficient,H=value%*%coefficient,g=g)
}

.np_bias_mixed_product <- function(b, W, h) {
  if(nrow(W)>500L)stop("mixed row block exceeded",call.=FALSE)
  out <- matrix(0,nrow(W),nrow(b$x))
  for(ii in .np_bias_profile_tiles(nrow(b$x)))
    out[,ii] <- W%*%.np_bias_mixed_kernel(b,b$x[ii,,drop=FALSE],h,TRUE)
  out
}

.np_bias_mixed_p0_owner <- function(payload, progress=FALSE) {
  force(payload)
  function(ii) {
    n <- nrow(payload$x)
    if(!length(ii) || length(ii)>16L || anyNA(ii) || any(ii<1L | ii>n | ii!=floor(ii)))
      stop("invalid mixed P0 row request",call.=FALSE)
    H <- matrix(0,length(ii),n)
    for(b in payload$blocks) for(j in which(ii %in% b$indices))
      H[j,b$indices] <- .np_bias_mixed_pilot(b,payload$x[ii[j],,drop=FALSE],b$geo$b0)$H
    if(isTRUE(progress)) .np_plot_progress_notify(done=max(ii),total=n,unit="training rows")
    H
  }
}

.np_bias_mixed_student_prepare <- function(state) {
  pilot <- list(mixed=state$mixed,x=state$xcon)
  rmu <- drop(.np_bias_residual_apply(.np_bias_mixed_p0_owner(state$mixed,progress=TRUE),
    matrix(state$mu,ncol=1L))$residual)
  list(pilot=pilot,q=state$residual_row_variance,rmu=rmu,root.variance=sqrt(state$variance),
    smoother=list(x=state$xcon,groups=state$groups,h=state$variance_h))
}

# All rows at an anchor share the same P0 and residual operator; P1 is composed
# once for the whole bounded batch. Main rows come from the captured fit owner.
.np_bias_mixed_batch <- function(s, grid, axes, main) {
  if(nrow(grid)>25L || nrow(grid)*length(axes)>500L)
    stop("mixed correction batch exceeded",call.=FALSE)
  orders <- vapply(axes,function(j).np_bias_derivative_order(s,j),integer(1))
  units <- vapply(seq_along(axes),function(k)if(axes[k])s$scale[axes[k]]^orders[k] else 1,numeric(1))
  at <- grid
  for(j in seq_len(s$d))at[[s$numeric.cols[j]]] <- (at[[s$numeric.cols[j]]]-s$offset[j])/s$scale[j]
  C <- matrix(0,nrow(grid)*length(axes),nrow(s$x))
  for(b in s$mixed$blocks) {
    geo <- b$geo; Q <- b$Q; W <- T1 <- matrix(0,nrow(C),nrow(Q))
    for(i in seq_len(nrow(at))) {
      e <- at[i,,drop=FALSE]; rr <- (i-1L)*length(axes)+seq_along(axes)
      po <- .np_bias_mixed_pilot(b,e,geo$b0); A <- po$A; g <- po$g
      present <- b$pooled || identical(.np_bias_profile_keys(e,b$fc)[1L],b$keys[1L])
      v <- do.call(rbind,lapply(seq_along(axes),function(k)
        .np_bias_mixed_basis(b,e,axes[k],orders[k])))*as.numeric(present)
      right <- function(z)z-(z%*%Q)%*%A
      kg <- function(z)sweep(.np_bias_mixed_product(b,z,geo$bK),2L,g,"*")
      projected <- sum(kg(A)*t(Q)); tau <- sum(g)-projected
      eta <- sum(g)-(1-geo$rho)*projected
      if(!is.finite(eta) || eta<=0 || tau < -1e-10*sum(g))
        stop("M1 invalid residual trace",call.=FALSE)
      T <- do.call(rbind,lapply(seq_along(axes),function(k) {
        z <- drop(.np_bias_mixed_kernel(b,e,geo$bK,TRUE))*as.numeric(present)
        j <- axes[k]
        if(j)z <- z*.np_bias_gaussian_derivatives((b$xcon[,j]-e[[b$nc[j]]])/geo$bK[j],geo$bK[j],orders[k])[,orders[k]+1L]
        z*g/eta
      }))
      L <- main[rr,b$indices,drop=FALSE]*units
      w <- L+(v-L%*%Q)%*%A
      delta <- right(T-kg(w)/eta)
      for(step in seq_len(s$mixed$iterations)) {
        if(step>1L)delta <- delta-right(kg(delta)/eta)
        w <- w+delta
      }
      W[rr,] <- w
      T1[rr,] <- do.call(rbind,lapply(seq_along(axes),function(k)
        .np_bias_mixed_pilot(b,e,geo$b1,axes[k],orders[k])$H))*as.numeric(present)
      .np_plot_progress_notify()
    }
    WP <- matrix(0,nrow(W),nrow(b$x))
    for(ii in .np_bias_profile_tiles(nrow(b$x))) {
      H <- do.call(rbind,lapply(ii,function(i).np_bias_mixed_pilot(b,b$x[i,,drop=FALSE],geo$b1)$H))
      WP <- WP+W[,ii,drop=FALSE]%*%H
      .np_plot_progress_notify(done=max(ii),total=nrow(b$x),unit="pilot rows")
    }
    C[,b$indices] <- (W-WP+T1)/rep(units,nrow(grid))
  }
  C
}

.np_bias_mixed_frames <- function(state, grid, map) {
  lapply(seq_len(nrow(map)),function(i)
    .np_bias_frames(state,grid[map$row[i],,drop=FALSE],map$target[i]))
}

.np_bias_mixed_match <- function(at, e) {
  ii <- which(vapply(seq_len(nrow(at)),function(j)
    all(vapply(seq_len(ncol(at)),function(k)identical(at[[k]][j],e[[k]][1L]),logical(1))),logical(1)))
  if(length(ii)!=1L)stop("mixed capture endpoint mismatch",call.=FALSE)
  ii
}

# Expand categorical endpoints BEFORE the one distributed fit request. The
# display keeps original rows; inference retains the complete compact capture.
.np_bias_mixed_eval <- function(state, bws, xdat, ydat, exdat, gradients, gradient.order, targets) {
  map <- expand.grid(row=seq_len(nrow(exdat)),target=targets,stringsAsFactors=FALSE)
  frames <- .np_bias_mixed_frames(state,exdat,map)
  endpoints <- do.call(rbind,lapply(frames,function(z)do.call(rbind,z$e)))
  at <- unique(rbind(exdat,endpoints)); rownames(at) <- NULL
  fit <- npreg(bws=bws,txdat=xdat,tydat=ydat,exdat=at,se=TRUE,gradients=gradients,
    gradient.order=gradient.order,warn.glp.gradient=FALSE,
    .np.defer.empty.rows=TRUE,.np.fit.inference=TRUE)
  index <- vapply(seq_len(nrow(exdat)),function(i).np_bias_mixed_match(at,exdat[i,,drop=FALSE]),integer(1))
  display <- fit
  for(name in c("mean","merr")) display[[name]] <- fit[[name]][index]
  for(name in c("grad","gerr")) if(!is.null(fit[[name]]))display[[name]] <- fit[[name]][index,,drop=FALSE]
  # Empty-row metadata uses the same display-row map as the displayed values.
  empty <- attr(fit,".np.empty.rows",exact=TRUE)
  if(!is.null(empty))attr(display,".np.empty.rows") <- empty[index]
  display$mixed.capture <- fit
  display
}

.np_bias_mixed_capture <- function(state, fit, grid, targets) {
  captured <- fit[["mixed.capture",exact=TRUE]]
  if(is.null(captured))captured <- fit
  request <- captured[["fit.request",exact=TRUE]]; expected <- state$main.request
  if(!is.list(request) || !.np_plot_regression_capture_same_bws(request$bws,expected$bws) ||
      !identical(request$xdat,expected$xdat) || !identical(as.double(request$ydat),expected$ydat))
    stop("mixed adaptation requires its complete fitted plot request",call.=FALSE)
  at <- request$exdat
  ridge <- captured[["fit.inference",exact=TRUE]][[3L]]
  if(length(ridge)!=nrow(at) || any(!is.finite(ridge)) || any(ridge!=0))
    stop("regularized fixed-bandwidth main rows remain unqualified for bias correction",call.=FALSE)
  frames <- .np_bias_mixed_frames(state,grid,expand.grid(row=seq_len(nrow(grid)),target=targets,stringsAsFactors=FALSE))
  axes <- sort(unique(vapply(frames,`[[`,integer(1),"nu")))
  width <- ncol(state$Q); n <- nrow(state$x); m <- nrow(at)
  max.bytes <- 256*(as.double(n)*min(500L,m)+(as.double(n)+m)*(state$d+width+1)+as.double(width)^2)
  operators <- lapply(axes,function(j) {
    nu <- integer(state$d); if(j)nu[j] <- .np_bias_derivative_order(state,j)
    derivative <- if(j)nu else NULL
    owner <- .np_plot_regression_capture_operator(captured,request$bws,request$xdat,
      request$ydat,at,derivative,NULL,max.bytes)
    list(owner=owner$operator,derivative=derivative)
  })
  function(index,map) {
    frames <- .np_bias_mixed_frames(state,grid[index,,drop=FALSE],map)
    endpoints <- unique(do.call(rbind,lapply(frames,function(z)do.call(rbind,z$e))))
    # Endpoint sums are accumulated as they are produced; no full-grid n rows.
    C <- L <- matrix(0,nrow(map),n)
    for(jj in .np_bias_profile_tiles(nrow(endpoints),size=min(25L,500L %/% length(axes)))) {
      e <- endpoints[jj,,drop=FALSE]
      captured.index <- vapply(seq_len(nrow(e)),function(i).np_bias_mixed_match(at,e[i,,drop=FALSE]),integer(1))
      main <- lapply(operators,function(op)op$owner(at[captured.index,,drop=FALSE],op$derivative,captured.index))
      rows <- do.call(rbind,lapply(seq_len(nrow(e)),function(i)do.call(rbind,lapply(main,function(z)z[i,,drop=FALSE]))))
      corrected <- .np_bias_mixed_batch(state,e,axes,rows)
      for(i in seq_along(frames))for(k in seq_along(frames[[i]]$e)) {
        ix <- .np_bias_mixed_match(endpoints,frames[[i]]$e[[k]])
        if(!ix %in% jj)next
        r <- (match(ix,jj)-1L)*length(axes)+match(frames[[i]]$nu,axes)
        C[i,] <- C[i,]+frames[[i]]$sign[k]*corrected[r,]
        L[i,] <- L[i,]+frames[[i]]$sign[k]*rows[r,]
      }
    }
    factor <- sweep(C,2L,sqrt(state$variance),"*")
    if(any(!is.finite(factor)))stop("nonfinite mixed adapted influence rows",call.=FALSE)
    list(W=C,original=L,center=drop(C%*%state$y),original_center=drop(L%*%state$y),
      se=sqrt(rowSums(factor^2)),factor=factor,
      main.regularization=data.frame(row=map$row,target=map$target,ridge=0),
      adaptation=list(version=state$mixed$version,geometry=state$mixed$geometry,
        iterations=state$mixed$iterations,variance_bandwidth=state$variance_h,
        residual_policy="adapted P0 row-variance normalized"))
  }
}
