# Coordinatewise local-to-global correction and public integration seam.
# Model fitting, calibration and variance smoothing retain their existing owners.
# Only bounded target rows and thin training/basis factors are retained.
.np_bias_adaptation_basis <- function(x, at, degree, derivative = 0L, order = 1L) {
  nu <- integer(ncol(x)); if (derivative) nu[derivative] <- order
  W.lp(as.data.frame(x), exdat=as.data.frame(at), degree=degree,
    gradient.vec=if (derivative) nu else NULL, basis="glp", bernstein.basis=TRUE)
}

.np_bias_adaptation_kernel <- function(x, at, bandwidth) {
  if (nrow(at)>16L) stop("internal adaptation kernel tile exceeded",call.=FALSE)
  out <- matrix(0,nrow(x),nrow(at))
  for (j in seq_along(bandwidth))
    out <- out+outer(x[,j],at[,j],"-")^2/bandwidth[j]^2
  exp(-out/2)
}

.np_bias_adaptation_product <- function(rows, x, bandwidth) {
  if (nrow(rows)>500L) stop("internal adaptation row block exceeded",call.=FALSE)
  out <- matrix(0,nrow(rows),nrow(x))
  for (first in seq.int(1L,nrow(x),by=16L)) {
    ii <- seq.int(first,min(nrow(x),first+15L))
    out[,ii] <- rows %*% .np_bias_adaptation_kernel(x,x[ii,,drop=FALSE],bandwidth)
    if ((first-1L) %% 256L == 0L) .np_plot_progress_notify()
  }
  out
}

.np_bias_adaptation_geometry <- function(p, h, b0, b1, bK, q) {
  if (length(p)!=length(q) || any(q<p) || any(!is.finite(h)) || any(h<=0))
    stop("invalid fixed-model adaptation geometry",call.=FALSE)
  r <- pmax(1,h/b0); u <- r^2-1; terms <- npBuildLpTerms(q)
  penalty <- numeric(nrow(terms))
  for (j in seq_along(q)) {
    excess <- pmax(terms[,j]-p[j],0); ii <- which(excess>0)
    penalty[ii] <- penalty[ii]+u[j]*excess[ii]^2
  }
  total <- rowSums(terms)
  for (k in unique(total[total>0])) {
    ii <- which(total==k)
    excess <- pmax(k-pmax(max(p),apply(terms[ii,,drop=FALSE],1,max)),0)
    jj <- which(excess>0)
    if (length(jj)) penalty[ii[jj]] <- penalty[ii[jj]]+min(u[q>=k])*excess[jj]^2
  }
  if (any(!is.finite(penalty))) stop("nonfinite adaptation penalty",call.=FALSE)
  list(q=q,r=r,lambda=penalty,rho=prod(1-r^-2),b0=b0*r,b1=b1*r,bK=bK*r,
    original.q=q,global=is.infinite(h),inactive=all(r==1))
}

# Q and its reference geometry are prepared once, not once per pilot row.
.np_bias_adaptation_pilot <- function(state, at, bandwidth, value = NULL,
                                    derivatives = integer(), orders = NULL) {
  .np_bias_transition_pilot(state, at, bandwidth, value, derivatives, orders)
}

.np_bias_adaptation_prepare <- function(spec) {
  if (length(spec$factor.cols) || spec$bws[["type"]]!="fixed" ||
      spec$bws[["ckertype"]]!="gaussian" || spec$bws[["ckerorder"]]!=2L ||
      spec$bws[["ckerbound"]]!="none" || !isTRUE(spec$bernstein))
    stop("private adaptation seam requires continuous fixed Gaussian GLP pilots with a stable basis",call.=FALSE)
  reg <- npValidatedConditionalRegSpec(spec$bws,where="private adaptation",ncon.field="ncon")
  if (reg[["basis.engine"]]!="glp") stop("private adaptation requires GLP",call.=FALSE)
  geo <- .np_bias_adaptation_geometry(reg[["degree.engine"]],
    unname(spec$bws[["bandwidth"]][["x"]]),spec$b,spec$p1$b,spec$h,spec$q)
  state <- list(spec=spec,x=as.matrix(spec$xcon),Q=spec$Q,geometry=geo,
    version="B4-nested-projection-1",iterations=spec$settings[["iterations"]])
  if (geo$inactive) {
    state$residual <- spec$residual; state$residual_raw <- spec$residual_raw
    state$residual_row_variance <- spec$residual_row_variance; state$variance <- spec$variance
    return(state)
  }
  state$geometry$transition <- .np_bias_transition_geometry(state$Q,geo$lambda,state$x)
  n <- nrow(state$x); raw <- rowvar <- numeric(n)
  for (i in seq_len(n)) {
    H <- .np_bias_adaptation_pilot(state,state$x[i,,drop=FALSE],geo$b0,
      state$Q[i,,drop=FALSE])$H
    R <- -drop(H); R[i] <- R[i]+1
    raw[i] <- sum(R*spec$y); rowvar[i] <- sum(R^2)
    if (i %% 16L == 0L || i == n)
      .np_plot_progress_notify(done=i,total=n,unit="training rows")
  }
  if (any(!is.finite(rowvar)) || any(rowvar<=0))
    stop("invalid adapted P0 residual row",call.=FALSE)
  state$residual_raw <- raw; state$residual_row_variance <- rowvar
  state$residual <- raw/sqrt(rowvar)
  state$variance <- .np_bias_variance(spec,state$residual)
  state
}

# Apply P1 once to every row in this bounded batch. No training-square pilot.
.np_bias_adaptation_compose <- function(rows, state) {
  x <- state$x; out <- matrix(0,nrow(rows),nrow(x))
  for (first in seq.int(1L,nrow(x),by=16L)) {
    ii <- seq.int(first,min(nrow(x),first+15L))
    H <- do.call(rbind,lapply(ii,function(i)
      .np_bias_adaptation_pilot(state,x[i,,drop=FALSE],state$geometry$b1,
        state$Q[i,,drop=FALSE])$H))
    out <- out+rows[,ii,drop=FALSE]%*%H
    .np_plot_progress_notify(done=max(ii),total=nrow(x),unit="pilot rows")
  }
  out
}

.np_bias_adaptation_rows <- function(state, grid, map, main.rows) {
  spec <- state$spec; x <- state$x; geo <- state$geometry; Q <- state$Q
  if (nrow(map)>500L || nrow(grid)>25L || nrow(main.rows)!=nrow(map) ||
      ncol(main.rows)!=nrow(x) || any(!is.finite(main.rows)))
    stop("invalid bounded adaptation target rows",call.=FALSE)
  if (geo$inactive) return(.np_bias_regression_rows(spec,grid,map))
  at <- sweep(sweep(as.matrix(grid),2L,spec$offset,"-"),2L,spec$scale,"/")
  nu <- vapply(map$target,function(t) if(t=="mean") 0L else {
    j <- suppressWarnings(as.integer(substring(t,2L)))
    if (!startsWith(t,"d") || is.na(j) || j<1L || j>spec$d)
      stop("invalid adaptation derivative target",call.=FALSE)
    j
  },integer(1))
  units <- vapply(nu,function(j)if(j)spec$scale[j]^.np_bias_derivative_order(spec,j) else 1,numeric(1))
  L <- main.rows*units; W <- T1 <- matrix(0,nrow(map),nrow(x))
  eta <- tau <- numeric(nrow(grid))
  for (anchor in unique(map$row)) {
    rr <- which(map$row==anchor); t <- at[anchor,,drop=FALSE]
    po <- .np_bias_adaptation_pilot(state,t,geo$b0); A <- po$A; g <- po$g
    v <- do.call(rbind,lapply(nu[rr],function(j).np_bias_adaptation_basis(x,t,geo$q,j,.np_bias_derivative_order(spec,j))))
    resid <- function(rows) rows-(rows%*%Q)%*%A
    kg <- function(rows) sweep(.np_bias_adaptation_product(rows,x,geo$bK),2L,g,"*")
    projected <- sum(diag(kg(A)%*%Q))
    tau[anchor] <- sum(g)-projected
    eta[anchor] <- sum(g)-(1-geo$rho)*projected
    normalizer <- eta[anchor]
    if (!is.finite(normalizer) || normalizer<=0 || tau[anchor] < -1e-10*sum(g))
      stop("invalid declared step normalization",call.=FALSE)
    T <- do.call(rbind,lapply(nu[rr],function(j) {
      v <- drop(.np_bias_adaptation_kernel(x,t,geo$bK))
      if(j) {
        order <- .np_bias_derivative_order(spec,j)
        if(order==1L)v <- v*(x[,j]-t[1,j])/geo$bK[j]^2
        else v <- v*.np_bias_gaussian_derivatives((x[,j]-t[1,j])/geo$bK[j],geo$bK[j],order)[,order+1L]
      }
      v*g/normalizer
    }))
    w <- L[rr,,drop=FALSE]+(v-L[rr,,drop=FALSE]%*%Q)%*%A
    delta <- resid(T-kg(w)/normalizer)
    for (step in seq_len(state$iterations)) {
      if(step>1L)delta <- delta-resid(kg(delta)/normalizer)
      w <- w+delta
    }
    W[rr,] <- w
    derivative.targets <- nu[rr][nu[rr]>0L]
    pilot <- .np_bias_adaptation_pilot(state,t,geo$b1,
      derivatives=derivative.targets,orders=spec$gradient.order)
    derivative.row <- 0L
    for (k in seq_along(rr)) {
      j <- nu[rr[k]]
      if(j) {
        derivative.row <- derivative.row+1L
        T1[rr[k],] <- pilot$D[derivative.row,]
      } else T1[rr[k],] <- pilot$H
    }
  }
  C <- (W-.np_bias_adaptation_compose(W,state)+T1)/units
  factor <- sweep(C,2L,sqrt(state$variance),"*")
  if (any(!is.finite(factor))) stop("nonfinite adapted influence rows",call.=FALSE)
  geo$transition <- NULL
  list(W=C,original=main.rows,center=drop(C%*%spec$y),
    original_center=drop(main.rows%*%spec$y),se=sqrt(rowSums(factor^2)),factor=factor,
    adaptation=list(version=state$version,geometry=geo,iterations=state$iterations,
      tau=tau,eta=eta,training_p1_rows=nrow(x),variance_bandwidth=spec$variance_h,
      residual_policy="adapted P0 row-variance normalized"))
}

# One activation policy for public callers. Old controls and literal overrides
# keep their original construction unless adaptation is explicitly requested.
.np_bias_adaptation_policy <- function(spec) {
  mode <- spec$settings$adaptation_request
  requested <- spec$settings$requested
  roles <- c("p0", "p1", "operator")
  defaults <- c(p0=4,p1=1.5,operator=4/3)
  override <- any(vapply(roles,function(role) {
    key <- paste0(role,".scale")
    value <- requested$bandwidth.controls[[key]]
    key %in% requested$supplied || is.null(value) || !identical(as.double(value),as.double(defaults[[role]])) ||
      !is.null(names(value)) || length(value)!=1L
  },logical(1)))
  eligible <- !length(spec$factor.cols) && spec$bws[["type"]]=="fixed" &&
    spec$bws[["ckertype"]]=="gaussian" && spec$bws[["ckerorder"]]==2L &&
    spec$bws[["ckerbound"]]=="none" && isTRUE(spec$bernstein)
  reason <- if(mode=="none") "disabled" else if(!eligible) "incumbent model class" else
    if(mode=="auto" && override) "explicit auxiliary geometry override" else "eligible"
  if(mode=="local-to-global" && !eligible)
    stop("local-to-global adaptation requires continuous-only fixed Gaussian GLP with stable auxiliary basis",call.=FALSE)
  active <- FALSE
  if(reason=="eligible") {
    active <- any(unname(spec$bws[["bandwidth"]][["x"]])>spec$b)
    if(!active) reason <- "local geometry unchanged"
  }
  list(requested=mode,eligible=eligible,active=active,reason=reason,
    version="B4-nested-projection-1")
}

# Capture is tied to the complete panel request, not a newly fitted row batch.
.np_bias_adaptation_capture <- function(state, fit, grid, targets) {
  request <- fit[["fit.request",exact=TRUE]]
  expected <- state$main.request
  if (!is.list(request) || !.np_plot_regression_capture_same_bws(request$bws,expected$bws) ||
      !identical(request$xdat,expected$xdat) || !identical(as.double(request$ydat),as.double(expected$ydat)) ||
      !identical(unname(as.matrix(request$exdat)),unname(as.matrix(grid))))
    stop("adapted bias correction requires its complete fitted plot request",call.=FALSE)
  ridge <- fit[["fit.inference",exact=TRUE]][[3L]]
  if (length(ridge)!=nrow(grid) || any(!is.finite(ridge)) || any(ridge!=0))
    stop("regularized fixed-bandwidth main rows remain unqualified for bias correction",call.=FALSE)
  width <- ncol(state$Q); n <- nrow(state$x); m <- nrow(grid)
  max.bytes <- 256*(as.double(n)*min(500L,m)+(as.double(n)+m)*(state$d+width+1)+as.double(width)^2)
  operators <- setNames(lapply(targets,function(target) {
    nu <- integer(state$d)
    if(target!="mean") {
      j <- suppressWarnings(as.integer(substring(target,2L)))
      if(!startsWith(target,"d") || is.na(j) || !j %in% seq_len(state$d))
        stop("invalid adaptation derivative target",call.=FALSE)
      nu[j] <- .np_bias_derivative_order(state,j)
    }
    derivative <- if(any(nu)) nu else NULL
    owner <- .np_plot_regression_capture_operator(fit,request$bws,request$xdat,
      request$ydat,request$exdat,derivative,NULL,max.bytes)
    list(owner=owner$operator,derivative=derivative)
  }),targets)
  function(index,map) {
    rows <- matrix(0,nrow(map),nrow(state$x))
    for(target in unique(map$target)) {
      rr <- which(map$target==target); ii <- index[map$row[rr]]; op <- operators[[target]]
      rows[rr,] <- op$owner(request$exdat[ii,,drop=FALSE],op$derivative,ii)
    }
    value <- .np_bias_adaptation_rows(state$adaptation,grid[index,,drop=FALSE],map,rows)
    value$main.regularization <- data.frame(row=map$row,target=map$target,frame=1L,ridge=ridge[index[map$row]])
    value
  }
}
