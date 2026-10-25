# Partially pooled auxiliary for newly admitted categorical layouts.
# Profile-specific main terms and shared higher curvature define one explicit
# polynomial space. No failure-driven fallback or training-square operator.
.np_bias_profile_layout <- function(bws, x) {
  fc <- which(!bws[["icon"]])
  length(fc) > 0L && (any(vapply(x[fc], nlevels, integer(1)) != 2L) ||
    sum(bws[["iuno"]]) > 1L || sum(bws[["iord"]]) > 1L)
}

.np_bias_profile_keys <- function(x, fc) {
  if (length(fc)) do.call(paste, c(lapply(x[fc], as.integer), sep=":"))
  else rep("all", nrow(x))
}

.np_bias_profile_tiles <- function(n, size=16L) {
  split(seq_len(n), (seq_len(n)-1L) %/% size)
}

.np_bias_profile_basis <- function(a, at, derivative=NULL) {
  Q <- W.lp(a$xcon, exdat=at[a$nc], degree=a$degree, gradient.vec=derivative, basis="glp",
    bernstein.basis=a$bernstein)
  M <- W.lp(a$xcon, exdat=at[a$nc], degree=a$main.degree, gradient.vec=derivative, basis="glp",
    bernstein.basis=a$bernstein)
  key <- .np_bias_profile_keys(at, a$fc)
  if (any(!key %in% a$levels))
    stop("pooled bias correction requires an observed evaluation profile", call.=FALSE)
  for (g in a$levels[-1L]) Q <- cbind(Q, M * as.numeric(key==g))
  Q
}

.np_bias_profile_solve <- function(Q, g, role="pooled auxiliary") {
  if (any(!is.finite(g)) || any(g<=0) || ncol(Q)>=nrow(Q))
    stop(role, " requires positive weights and residual support", call.=FALSE)
  Z <- Q * sqrt(g)
  norms <- sqrt(colSums(Z^2))
  if (any(!is.finite(norms)) || any(norms<=0))
    stop(role, " has an unsupported auxiliary column", call.=FALSE)
  z <- qr(sweep(Z, 2L, norms, "/"), tol=1e-11, LAPACK=FALSE)
  if (z$rank != ncol(Q))
    stop(role, " is numerically unidentified (rank ", z$rank, "/", ncol(Q),
      "); inspect local profile support and pilot bandwidth/degree controls", call.=FALSE)
  U <- qr.Q(z)
  A <- backsolve(qr.R(z), t(U)*rep(sqrt(g), each=ncol(Q)))
  coefficient <- matrix(0, ncol(Q), nrow(Q))
  coefficient[z$pivot,] <- A
  list(coefficient=sweep(coefficient, 1L, norms, "/"), thin=U*sqrt(g))
}

.np_bias_profile_kernel <- function(a, at, h, category=TRUE) {
  if (nrow(at)>16L) stop("pooled auxiliary kernel tile exceeded", call.=FALSE)
  K <- matrix(1, nrow(a$x), nrow(at))
  for (j in seq_along(h))
    K <- K * exp(-.5*(outer(a$xcon[[j]], at[[a$nc[j]]], "-")/h[j])^2)
  if (category) K <- K * outer(a$keys, .np_bias_profile_keys(at,a$fc), "==")
  K
}

.np_bias_profile_pilot_rows <- function(a, at, h, role) {
  E <- .np_bias_profile_basis(a, at)
  K <- .np_bias_profile_kernel(a, at, h, FALSE)
  out <- matrix(0, nrow(at), nrow(a$x))
  for (j in seq_len(nrow(at)))
    out[j,] <- drop(E[j,,drop=FALSE] %*%
      .np_bias_profile_solve(a$Q, K[,j], role)$coefficient)
  out
}

# The worker payload contains data/basis geometry, never a serialized closure
# capturing a package namespace, fitted response, or the entire inference state.
.np_bias_profile_pilot <- function(a) {
  a[c("x","xcon","nc","fc","keys","levels","degree","main.degree",
    "bernstein","Q","h")]
}

.np_bias_profile_prepare <- function(bws, x, y, control, reg, orders=NULL) {
  if (!bws[["type"]] %in% c("fixed", "adaptive_nn", "generalized_nn"))
    stop("unsupported pooled main bandwidth type", call.=FALSE)
  nc <- which(bws[["icon"]]); fc <- which(!bws[["icon"]])
  keys <- .np_bias_profile_keys(x,fc)
  groups <- split(seq_len(nrow(x)),keys)
  settings <- .np_bias_resolve(control,x[nc],nrow(x),reg[["degree.engine"]])
  if (settings$adaptation_request=="local-to-global")
    stop("local-to-global adaptation remains unqualified for pooled mixed-data auxiliaries", call.=FALSE)
  # Retain the frozen global-then-profile scale operation order. Raw overrides
  # are literal predictor-unit bandwidths and must not receive this adjustment.
  adjustment <- (nrow(x)/min(lengths(groups)))^settings$exponent
  for (role in names(settings$bandwidth))
    if (settings$source[[role]]=="profile_scale")
      settings$bandwidth[[role]] <- settings$bandwidth[[role]]*adjustment
  settings$profile_sizes <- lengths(groups)
  settings$n_min <- min(lengths(groups))
  settings$policy <- "pooled-curvature-v1"
  settings$effective_package_scale <- lapply(settings$bandwidth,
    function(h) h/(settings$spread*nrow(x)^(-settings$exponent)))
  degree <- unname(settings$degree_spec$effective)
  main.degree <- reg[["degree.engine"]]
  main.terms <- npLpCompleteTermStatus("glp",main.degree,nobs=min(lengths(groups)))
  width <- settings$degree_spec$nterms + (length(groups)-1)*main.terms$nterms
  if (!main.terms$capacity.ok || !main.terms$rank.ok || width>=nrow(x))
    stop("pooled auxiliary requires residual support in every profile main space and in the joint space", call.=FALSE)
  if (width>384L)
    stop("pooled auxiliary exceeds the qualified 384-column workspace; existing profile-specific routes are unchanged", call.=FALSE)
  offset <- vapply(x[nc],mean,numeric(1)); scale <- settings$spread
  xs <- x
  for (j in seq_along(nc)) xs[[nc[j]]] <- (x[[nc[j]]]-offset[j])/scale[j]
  a <- list(pooled=TRUE,bw=bws,x=xs,original.x=x,xcon=xs[nc],nc=nc,fc=fc,
    numeric.cols=nc,factor.cols=fc,d=length(nc),gradient.order=orders,keys=keys,levels=names(groups),groups=groups,
    main.degree=main.degree,degree=degree,bernstein=settings$bernstein.basis,
    offset=offset,scale=scale,y=as.double(y),settings=settings)
  a$h <- lapply(settings$bandwidth,function(h)unname(h/scale))
  a$Q <- .np_bias_profile_basis(a,xs)
  for (ii in groups) {
    M <- W.lp(a$xcon[ii,,drop=FALSE],degree=main.degree,basis="glp",
      bernstein.basis=a$bernstein)
    if (nrow(M)<=ncol(M) || qr(M)$rank!=ncol(M))
      stop("profile main space lacks identified residual support", call.=FALSE)
  }
  pilot <- .np_bias_profile_pilot(a)
  raw <- q <- mu <- numeric(nrow(x))
  for (ii in .np_bias_profile_tiles(nrow(x))) {
    H <- .np_bias_profile_pilot_rows(pilot,pilot$x[ii,,drop=FALSE],pilot$h$p0,"pooled P0")
    mu[ii] <- drop(H %*% y); raw[ii] <- y[ii]-mu[ii]
    H[cbind(seq_along(ii),ii)] <- H[cbind(seq_along(ii),ii)]-1
    q[ii] <- rowSums(H^2)
  }
  if (any(!is.finite(q)) || any(q<=0))
    stop("nonpositive pooled residual row variance",call.=FALSE)
  a$residual_raw <- raw; a$residual_row_variance <- q; a$mu <- mu
  a$residual <- raw/sqrt(q); a$variance_h <- a$h$variance
  smoother <- list(x=a$xcon,groups=a$groups,h=a$variance_h)
  S <- .np_bias_student_variance_owner(smoother)
  a$variance <- numeric(nrow(x))
  for (ii in .np_bias_profile_tiles(nrow(x)))
    a$variance[ii] <- drop(S(ii) %*% a$residual^2)
  if (any(!is.finite(a$variance)) || any(a$variance<=0))
    stop("nonpositive pooled smoothed variance",call.=FALSE)
  a$settings$auxiliary_space <- list(construction="shared curvature with profile-main deviations",
    columns=ncol(a$Q),shared.columns=settings$degree_spec$nterms,
    profile.main.columns=main.terms$nterms,profiles=names(groups),rank.tolerance=1e-11,
    workspace.columns=384L)
  if (!is.null(orders)) a$settings$gradient.order <- setNames(orders,names(a$xcon))
  a$settings$coordinate_offset <- offset; a$settings$coordinate_scale <- scale
  a$settings$training <- list(n=nrow(x),columns=names(x),row.names=rownames(x),
    rows.omit=bws[["rows.omit"]],source="retained plot training frame")
  a$settings$main <- list(bandwidth=bws[["bandwidth"]][["x"]],type=bws[["type"]],
    regtype=bws[["regtype"]],degree=bws[["degree"]],basis=bws[["basis"]],
    bernstein.basis=bws[["bernstein.basis"]],kernel=bws[["ckertype"]],order=bws[["ckerorder"]],
    scale.factor=bws[["sfactor"]][["x"]])
  a$settings$adaptation <- list(requested=settings$adaptation_request,active=FALSE,
    eligible=FALSE,reason="pooled mixed-data auxiliary")
  a
}

.np_bias_profile_student_prepare <- function(a) {
  pilot <- .np_bias_profile_pilot(a)
  rmu <- numeric(nrow(a$x))
  for (ii in .np_bias_profile_tiles(nrow(a$x)))
    rmu[ii] <- (a$mu[ii]-drop(.np_bias_profile_pilot_rows(pilot,
      pilot$x[ii,,drop=FALSE],pilot$h$p0,"pooled P0") %*% a$mu))/sqrt(a$residual_row_variance[ii])
  list(pilot=list(x=a$xcon,profile=pilot),q=a$residual_row_variance,rmu=rmu,
    root.variance=sqrt(a$variance),smoother=list(x=a$xcon,groups=a$groups,h=a$variance_h))
}

.np_bias_profile_product <- function(W, a, weight) {
  out <- matrix(0,nrow(W),ncol(W))
  for (ii in .np_bias_profile_tiles(nrow(a$x)))
    out[,ii] <- W %*% .np_bias_profile_kernel(a,a$x[ii,,drop=FALSE],a$h$operator)
  sweep(out,2L,weight,"*")
}

.np_bias_profile_component_rows <- function(a, grid, map, axes=integer(nrow(map))) {
  if (nrow(map)>25L || length(axes)!=nrow(map))
    stop("pooled mixed-data correction requires a bounded functional block",call.=FALSE)
  grid <- grid[map$row,,drop=FALSE]
  if (all(axes==0L)) {
    L <- npreghat(a$bw,txdat=a$original.x,exdat=grid,output="matrix",.np.require.finite=TRUE)
    ridge <- attr(L,"ridge.used",exact=TRUE)
  } else {
    L <- matrix(0,nrow(grid),nrow(a$x)); ridge <- numeric(nrow(grid))
    for (axis in unique(axes)) {
      ii <- which(axes==axis); nu <- integer(a$d)
      if (axis) nu[axis] <- .np_bias_derivative_order(a,axis)
      H <- npreghat(a$bw,txdat=a$original.x,exdat=grid[ii,,drop=FALSE],
        s=nu,output="matrix",.np.require.finite=TRUE)
      L[ii,] <- H; ridge[ii] <- attr(H,"ridge.used",exact=TRUE)
    }
  }
  if (length(ridge)!=nrow(grid) || any(!is.finite(ridge)) || any(ridge!=0))
    stop("pooled mixed-data correction requires unregularized main rows",call.=FALSE)
  attr(L,"ridge.used") <- NULL
  e <- grid
  for (j in seq_along(a$nc)) e[[a$nc[j]]] <- (e[[a$nc[j]]]-a$offset[j])/a$scale[j]
  E <- .np_bias_profile_basis(a,e)
  units <- rep(1,nrow(e))
  for (axis in unique(axes[axes>0L])) {
    ii <- which(axes==axis); nu <- integer(a$d)
    nu[axis] <- .np_bias_derivative_order(a,axis)
    units[ii] <- a$scale[axis]^nu[axis]
    E[ii,] <- .np_bias_profile_basis(a,e[ii,,drop=FALSE],nu)/units[ii]
  }
  W <- matrix(0,nrow(e),nrow(a$x))
  for (r in seq_len(nrow(e))) {
    g <- drop(.np_bias_profile_kernel(a,e[r,,drop=FALSE],a$h$p0,FALSE))
    p <- .np_bias_profile_solve(a$Q,g,"pooled P0 projection")
    projected <- 0
    for (ii in .np_bias_profile_tiles(nrow(a$x)))
      projected <- projected+sum(p$thin[ii,,drop=FALSE]*crossprod(
        .np_bias_profile_kernel(a,a$x[ii,,drop=FALSE],a$h$operator),p$thin))
    D <- sum(g)-projected
    if (!is.finite(D) || D<=1e-8*sum(g))
      stop("insufficient pooled residual trace",call.=FALSE)
    right <- function(Z) Z-(Z %*% a$Q) %*% p$coefficient
    w <- L[r,,drop=FALSE]+(E[r,,drop=FALSE]-L[r,,drop=FALSE] %*% a$Q) %*% p$coefficient
    target <- t(.np_bias_profile_kernel(a,e[r,,drop=FALSE],a$h$operator))*rep(g/D,each=1L)
    if (axes[r]) {
      axis <- axes[r]; order <- .np_bias_derivative_order(a,axis)
      target <- target * rep(.np_bias_gaussian_derivatives(
        (a$xcon[[axis]]-e[[a$nc[axis]]][r])/a$h$operator[axis],
        a$h$operator[axis],order)[,order+1L]/units[r],each=1L)
    }
    delta <- right(target-.np_bias_profile_product(w,a,g/D))
    for (step in seq_len(a$settings$iterations)) {
      if (step>1L) delta <- delta-right(.np_bias_profile_product(delta,a,g/D))
      w <- w+delta
    }
    W[r,] <- w
  }
  pilot <- .np_bias_profile_pilot(a)
  WP <- .np_bias_right_compose(W,function(ii)
    .np_bias_profile_pilot_rows(pilot,pilot$x[ii,,drop=FALSE],pilot$h$p1,"pooled P1"))
  T <- do.call(rbind,lapply(.np_bias_profile_tiles(nrow(e)),function(ii)
    .np_bias_profile_pilot_rows(pilot,e[ii,,drop=FALSE],pilot$h$p1,"pooled P1")))
  for (r in which(axes>0L)) {
    axis <- axes[r]; order <- .np_bias_derivative_order(a,axis)
    T[r,] <- .np_bias_profile_pilot_derivative(a,e[r,,drop=FALSE],axis,order)/units[r]
  }
  C <- W-WP+T
  F <- sweep(C,2L,sqrt(a$variance),"*"); se <- sqrt(rowSums(F^2))
  if (any(!is.finite(F)) || any(!is.finite(se)) || any(se<=0))
    stop("invalid pooled corrected influence rows",call.=FALSE)
  list(W=C,original=L,center=drop(C %*% a$y),original_center=drop(L %*% a$y),
    se=se,factor=F,main.regularization=data.frame(row=map$row,target="mean",frame=1L,ridge=ridge))
}


# Total derivative of the weighted joint-space pilot in standardized units.
.np_bias_profile_pilot_derivative <- function(a, at, axis, order) {
  value <- do.call(rbind,lapply(0:order,function(k) {
    nu <- integer(length(a$nc)); nu[axis] <- k
    .np_bias_profile_basis(a,at,if(k)nu else NULL)
  }))
  h <- a$h$p1
  weight <- drop(.np_bias_profile_kernel(a,at,h,FALSE))
  ratios <- .np_bias_gaussian_derivatives(
    (a$xcon[[axis]]-at[[a$nc[axis]]])/h[axis],h[axis],order)
  .np_bias_weighted_derivative(a$Q,weight,value,ratios,
    role="pooled P1 derivative",rank.tolerance=1e-11)
}

# Combine complete observation rows before covariance/studentization. In
# particular, categorical endpoints share their observation innovations.
.np_bias_profile_rows <- function(a, grid, map) {
  if (all(map$target=="mean"))
    return(.np_bias_profile_component_rows(a,grid,map))
  if (nrow(map)>500L) stop("pooled functional block exceeded",call.=FALSE)
  frames <- lapply(seq_len(nrow(map)),function(r)
    .np_bias_frames(a,grid[map$row[r],,drop=FALSE],map$target[r]))
  at <- do.call(rbind,lapply(frames,function(f)do.call(rbind,f$e)))
  owner <- rep(seq_along(frames),lengths(lapply(frames,`[[`,"e")))
  signs <- unlist(lapply(frames,`[[`,"sign"))
  axes <- rep(vapply(frames,`[[`,integer(1),"nu"),lengths(lapply(frames,`[[`,"e")))
  C <- L <- matrix(0,nrow(map),nrow(a$x)); regularization <- list()
  for (ii in .np_bias_profile_tiles(nrow(at),25L)) {
    value <- .np_bias_profile_component_rows(a,at[ii,,drop=FALSE],
      data.frame(row=seq_along(ii),target=map$target[owner[ii]]),axes[ii])
    for (j in seq_along(ii)) {
      r <- owner[ii[j]]
      C[r,] <- C[r,,drop=FALSE]+signs[ii[j]]*value$W[j,,drop=FALSE]
      L[r,] <- L[r,,drop=FALSE]+signs[ii[j]]*value$original[j,,drop=FALSE]
    }
    meta <- value$main.regularization
    meta$row <- map$row[owner[ii]]; meta$target <- map$target[owner[ii]]
    meta$frame <- sequence(lengths(lapply(frames,`[[`,"e")))[ii]
    regularization[[length(regularization)+1L]] <- meta
  }
  F <- sweep(C,2L,sqrt(a$variance),"*"); se <- sqrt(rowSums(F^2))
  if (any(!is.finite(F)) || any(!is.finite(se)))
    stop("invalid pooled corrected functional rows",call.=FALSE)
  list(W=C,original=L,center=drop(C%*%a$y),original_center=drop(L%*%a$y),
    se=se,factor=F,main.regularization=do.call(rbind,regularization))
}

# An automatic observed-profile inference grid may omit unused values, but its
# factor levels must remain intact: they define the fitted categorical kernel
# and the original reference/adjacent contrast endpoints.
.np_bias_observed_factor_values <- function(x) {
  if (!is.factor(x)) stop("observed categorical grid requires a factor",call.=FALSE)
  values <- levels(x)
  values[values %in% as.character(x[!is.na(x)])]
}
