# Regression adapter for the complete corrected observation operator.
# Training work is prepared once; each query batch has at most 25 anchors.
.np_bias_regbw <- function(source, x, bandwidth, degree = NULL, bernstein = NULL) {
  source[["scaling"]] <- FALSE
  if (!is.null(degree)) {
    source[["regtype"]] <- "lp"
    source[["basis"]] <- "glp"
    source[["degree"]] <- npValidateGlpDegree("lp", degree, ncol(x), argname="auxiliary degree")
    source[["bernstein.basis"]] <- bernstein
    source[["type"]] <- "fixed"
    source[["ckertype"]] <- "gaussian"
    source[["ckerorder"]] <- 2L
    source[["ckerbound"]] <- "none"
    source[["ckerlb"]] <- source[["ckerub"]] <- NULL
  }
  out <- .np_semihat_make_regbw_state(source, x, numeric(nrow(x)), bandwidth)
  out[["call"]] <- NULL
  out
}

.np_bias_regression_prepare <- function(bws, x, y, control = np_bias_control(), adapt = FALSE, gradient.order = NULL, law = NULL) {
  x <- toFrame(x)
  if (!is.numeric(y) || length(y) != nrow(x) || any(!is.finite(y)) || anyNA(x))
    stop("bias correction requires finite retained training data", call. = FALSE)
  reg <- npValidatedConditionalRegSpec(bws, where = "bias correction", ncon.field = "ncon")
  nc <- which(bws[["icon"]]); fc <- which(!bws[["icon"]]); d <- length(nc)
  if (!d) return(.np_bias_categorical_prepare(bws,x,y,control,reg))
  legacy.hint <- if (bws[["type"]]=="fixed" && d>0L)
    "; explicit bootstrap bias.control=np_bias_control(construction='legacy') retains the previous construction for supported bootstrap methods" else ""
  if (!d || bws[["ckertype"]] != "gaussian" || bws[["ckerorder"]] != 2L ||
      bws[["ckerbound"]] != "none" || reg[["basis.engine"]] != "glp")
    stop(paste0("iterated regression bias correction currently supports unbounded second-order Gaussian kernels, at least one continuous predictor and GLP main models",legacy.hint), call. = FALSE)
  if (length(fc) && !all(vapply(x[fc],is.factor,logical(1))))
    stop("bias correction requires factor-valued categorical predictors",call.=FALSE)
  orders <- .np_bias_derivative_request(reg,d,gradient.order)
  if (.np_bias_profile_layout(bws,x)) {
    return(.np_bias_profile_prepare(bws,x,y,control,reg,orders,adapt,law))
  }
  keys <- if (length(fc)) do.call(paste,c(lapply(x[fc],as.integer),sep=":")) else rep("all",nrow(x))
  groups <- split(seq_len(nrow(x)),keys)
  if (length(groups) != 2^length(fc)) stop("bias correction requires every training profile",call.=FALSE)
  settings <- .np_bias_resolve(control,x[nc],lengths(groups),reg[["degree.engine"]])
  q <- unname(settings[["degree_spec"]][["effective"]])
  for (ii in groups) {
    distinct <- vapply(x[ii,nc,drop=FALSE],function(z)length(unique(z)),integer(1))
    if (any(distinct <= q))
      stop("auxiliary polynomial degree requires more distinct values in each training profile; use bias.control=np_bias_control(degree=...) with an identified space containing the fixed main model",call.=FALSE)
  }
  offset <- vapply(x[nc],mean,numeric(1)); scale <- settings[["spread"]]
  standardized <- x
  for (j in seq_len(d)) standardized[[nc[j]]] <- (x[[nc[j]]]-offset[j])/scale[j]
  h <- unname(bws[["bandwidth"]][["x"]])
  if (bws[["type"]] == "fixed") h[nc] <- h[nc]/scale
  main <- .np_bias_regbw(bws,standardized,h)
  spec <- list(x=standardized,xcon=standardized[nc],numeric.cols=nc,factor.cols=fc,
    d=d,q=q,bernstein=settings[["bernstein.basis"]],keys=keys,groups=groups,bws=main,
    b=unname(settings[["bandwidth"]][["p0"]]/scale),
    h=unname(settings[["bandwidth"]][["operator"]]/scale),
    variance_h=unname(settings[["bandwidth"]][["variance"]]/scale),
    offset=offset,scale=scale,y=as.double(y),settings=settings)
  spec$pilot.bws <- .np_bias_regbw(bws,spec$xcon,spec$b,q,spec$bernstein)
  spec$Q <- W.lp(spec$xcon,degree=q,basis="glp",bernstein.basis=spec$bernstein)
  spec$p1 <- list(b=unname(settings[["bandwidth"]][["p1"]]/scale))
  spec$p1$pilot.bws <- .np_bias_regbw(bws,spec$xcon,spec$p1$b,q,spec$bernstein)
  spec$settings$coordinate_offset <- offset
  spec$settings$coordinate_scale <- scale
  spec$settings$training <- list(n=nrow(x),columns=names(x),row.names=rownames(x),
    rows.omit=bws[["rows.omit"]],source="retained plot training frame")
  spec$settings$main <- list(bandwidth=bws[["bandwidth"]][["x"]],type=bws[["type"]],
    regtype=bws[["regtype"]],degree=bws[["degree"]],basis=bws[["basis"]],
    bernstein.basis=bws[["bernstein.basis"]],kernel=bws[["ckertype"]],order=bws[["ckerorder"]],
    scale.factor=bws[["sfactor"]][["x"]])
  if (!is.null(orders)) {
    spec$gradient.order <- orders
    spec$settings$gradient.order <- setNames(orders,names(spec$xcon))
  }
  policy <- if (isTRUE(adapt)) {
    if(length(fc)) .np_bias_mixed_policy(spec,bws,law) else .np_bias_adaptation_policy(spec)
  } else NULL
  if(length(fc) && !is.null(policy) && policy$active)
    return(.np_bias_mixed_prepare(spec,bws,x,policy))
  if (!is.null(policy) && policy$active) {
    adapted <- .np_bias_adaptation_prepare(spec)
    for (name in c("residual", "residual_raw", "residual_row_variance", "variance"))
      spec[[name]] <- adapted[[name]]
    spec$adaptation <- adapted
    spec$main.request <- list(bws=bws,xdat=x,ydat=as.double(y))
    spec$settings$policy <- "rot-profile-v1+continuous-nested-projection-v1"
    spec$settings$reference_bandwidth <- spec$settings$bandwidth
    spec$settings$reference_package_scale <- spec$settings$effective_package_scale
    effective <- list(p0=adapted$geometry$b0,p1=adapted$geometry$b1,
      operator=adapted$geometry$bK,variance=spec$variance_h)
    for (role in names(effective)) {
      spec$settings$bandwidth[[role]] <- setNames(effective[[role]]*scale,names(spec$xcon))
      spec$settings$effective_package_scale[[role]] <-
        spec$settings$bandwidth[[role]]/(scale*nrow(x)^(-spec$settings$exponent))
    }
    policy$geometry <- adapted$geometry
    policy$geometry$transition <- NULL
    policy$residual_policy <- "adapted P0 row-variance normalized"
  } else {
    training <- .np_bias_p0_residuals(spec,matrix(y,ncol=1L))
    spec$residual <- drop(training$residual)
    spec$residual_raw <- drop(training$raw)
    spec$residual_row_variance <- training$row_variance
    spec$variance <- .np_bias_variance(spec,spec$residual)
  }
  if (!is.null(policy)) spec$settings$adaptation <- policy
  spec
}

.np_bias_variance <- function(spec, residual, block = 16L) {
  bws <- .np_bias_regbw(spec$pilot.bws,spec$xcon,spec$variance_h,rep.int(1L,spec$d),FALSE)
  variance <- numeric(nrow(spec$x))
  for (ii in spec$groups) for (first in seq.int(1L,length(ii),by=block)) {
    jj <- seq.int(first,min(length(ii),first+block-1L))
    K <- .np_kernel_weights_direct(bws,txdat=spec$xcon[ii,,drop=FALSE],
      exdat=spec$xcon[ii[jj],,drop=FALSE],bandwidth.divide=TRUE)
    mass <- colSums(K)
    if (any(!is.finite(K)) || any(K<0) || any(mass<=0)) stop("invalid bias variance smoother",call.=FALSE)
    H <- t(sweep(K,2L,mass,"/"))
    variance[ii[jj]] <- drop(H %*% residual[ii]^2)
  }
  if (any(!is.finite(variance)) || any(variance<=0)) stop("nonpositive corrected variance",call.=FALSE)
  variance
}

.np_bias_frames <- function(spec,e,target) {
  if (startsWith(target,"c")) {
    j <- as.integer(substring(target,2L))
    if (!j %in% spec$factor.cols) stop("invalid categorical bias target",call.=FALSE)
    f <- npCategoricalFirstDifferenceFrames(e,j,"bias correction")
    list(e=list(f$upper,f$lower),sign=c(1,-1),nu=0L)
  } else {
    nu <- if (target=="mean") 0L else as.integer(substring(target,2L))
    if (is.na(nu) || !nu %in% 0:spec$d) stop("invalid derivative bias target",call.=FALSE)
    list(e=list(e),sign=1,nu=nu,order=.np_bias_derivative_order(spec,nu))
  }
}

.np_bias_projector_blocks <- function(spec,anchor,reference=FALSE) {
  if (reference) stop("runtime bias operator has no diagnostic fallback",call.=FALSE)
  .np_bias_projector(spec,anchor)$blocks
}

.np_bias_main_rows <- function(spec,evals,nu) {
  s <- integer(spec$d); if(nu)s[nu] <- .np_bias_derivative_order(spec,nu)
  reg <- npValidatedConditionalRegSpec(spec$bws,where="bias main operator",ncon.field="ncon")
  if (spec$bws[["type"]]=="generalized_nn" && any(reg[["degree.engine"]]>1L)) {
    # Retain the qualified prepared-adjoint GNN owner used by the checkpoint.
    K <- .np_kernel_weights_direct(spec$bws,txdat=spec$x,exdat=evals,bandwidth.divide=TRUE)
    P <- W.lp(spec$xcon,degree=reg[["degree.engine"]],basis=reg[["basis.engine"]],bernstein.basis=reg[["bernstein.basis.engine"]])
    V <- W.lp(spec$xcon,exdat=evals[spec$numeric.cols],degree=reg[["degree.engine"]],
      gradient.vec=if(nu)s else NULL,basis=reg[["basis.engine"]],bernstein.basis=reg[["bernstein.basis.engine"]])
    H <- matrix(0,nrow(evals),nrow(spec$x));ridge <- numeric(nrow(evals))
    for(j in seq_len(nrow(evals))) {
      weighted <- P*K[,j]
      solved <- .Call("C_np_lp_adjoint_prepared_ridge",crossprod(P,weighted),as.double(V[j,]),
        as.integer(nrow(P)),as.integer(min(ncol(P),sum(K[,j]!=0))),PACKAGE="npRmpi")
      ridge[j] <- solved$ridge
      H[j,] <- drop(weighted %*% solved$v)
    }
  } else {
    H <- npreghat(spec$bws,txdat=spec$x,exdat=evals,s=s,output="matrix",.np.require.finite=TRUE)
    ridge <- attr(H,"ridge.used",exact=TRUE)
  }
  if (is.null(ridge) || any(!is.finite(ridge)) || any(ridge<0) || any(!is.finite(H)))
    stop("bias correction requires finite canonical main observation rows",call.=FALSE)
  admitted <- spec$bws[["type"]]=="generalized_nn" &&
    any(reg[["degree.engine"]]>1L) && isTRUE(reg[["bernstein.basis.engine"]])
  if (any(ridge!=0) && !admitted)
    stop("regularized main rows are qualified only for generalized-NN GLP models with a stored Bernstein basis; the fitted model has not been changed",call.=FALSE)
  attr(H,"ridge.used") <- ridge
  H
}

.np_bias_original <- function(spec,grid,map,block=16L) {
  out <- matrix(0,nrow(map),nrow(spec$x))
  regularization <- list()
  for(target in unique(map$target)) {
    rr <- which(map$target==target)
    for(first in seq.int(1L,length(rr),by=block)) {
      jj <- rr[seq.int(first,min(length(rr),first+block-1L))]
      frames <- lapply(jj,function(j).np_bias_frames(spec,grid[map$row[j],,drop=FALSE],target))
      for(k in seq_along(frames[[1L]]$sign)) {
        evals <- do.call(rbind,lapply(frames,function(f)f$e[[k]]))
        H <- .np_bias_main_rows(spec,evals,frames[[1L]]$nu)
        regularization[[length(regularization)+1L]] <- data.frame(
          row=map$row[jj],target=target,frame=k,ridge=attr(H,"ridge.used",exact=TRUE),
          stringsAsFactors=FALSE)
        out[jj,] <- out[jj,,drop=FALSE]+frames[[1L]]$sign[k]*H
      }
    }
  }
  attr(out,"main.regularization") <- do.call(rbind,regularization)
  out
}

.np_bias_polynomial <- function(spec,grid,map,L,projectors) {
  W <- L
  for(a in projectors)for(b in a$blocks) {
    rr <- a$rows; ii <- b$indices; key <- spec$keys[ii[1L]]
    DQ <- t(vapply(rr,function(r) {
      f <- .np_bias_frames(spec,grid[map$row[r],,drop=FALSE],map$target[r]);v <- numeric(ncol(spec$Q))
      for(k in seq_along(f$e))if(.np_bias_key(spec,f$e[[k]])==key) {
        s <- integer(spec$d);if(f$nu)s[f$nu] <- f$order
        v <- v+f$sign[k]*.np_bias_projection(spec,f$e[[k]],spec$q,s)
      }
      v
    },numeric(ncol(spec$Q))))
    W[rr,ii] <- L[rr,ii,drop=FALSE]+(DQ-L[rr,ii,drop=FALSE]%*%b$Q)%*%b$coefficient
  }
  W
}

.np_bias_regression_rows <- function(spec,grid,map) {
  if (isTRUE(spec$categorical.only)) return(.np_bias_categorical_rows(spec,grid,map))
  if (isTRUE(spec$pooled)) return(.np_bias_profile_rows(spec,grid,map))
  if (nrow(map)>500L || nrow(unique(grid[spec$numeric.cols]))>25L)
    stop("internal bias target block exceeded",call.=FALSE)
  for(j in seq_len(spec$d))grid[[spec$numeric.cols[j]]] <- (grid[[spec$numeric.cols[j]]]-spec$offset[j])/spec$scale[j]
  L <- .np_bias_original(spec,grid,map)
  regularization <- attr(L,"main.regularization",exact=TRUE)
  attr(L,"main.regularization") <- NULL
  op <- .np_bias_stable_prepare(spec,grid,map)
  W <- .np_bias_polynomial(spec,grid,map,L,op$projectors)
  target <- .np_bias_stable_target(spec,grid,map,op)
  delta <- .np_bias_right_residual(target-.np_bias_stable_product(W,op),op$projectors)
  for(step in seq_len(spec$settings[["iterations"]])) {
    if(step>1L)delta <- delta-.np_bias_right_residual(.np_bias_stable_product(delta,op),op$projectors)
    W <- W+delta
  }
  pilot <- spec;pilot$b <- spec$p1$b;pilot$pilot.bws <- spec$p1$pilot.bws
  WP <- .np_bias_right_compose(W,.np_bias_pilot_owner(pilot))
  # Category contrasts use differences of complete pilot rows, including covariance.
  T <- matrix(0,nrow(map),nrow(spec$x))
  for(r in seq_len(nrow(map))) {
    f <- .np_bias_frames(spec,grid[map$row[r],,drop=FALSE],map$target[r])
    for(k in seq_along(f$e)) {
      p <- .np_bias_pilot_rows(pilot,f$e[[k]],derivatives=f$nu>0L,orders=spec$gradient.order)
      h <- if(f$nu)p$derivatives[[f$nu]] else p$mean
      T[r,p$indices] <- T[r,p$indices]+f$sign[k]*h
    }
  }
  C <- W-(WP-T)
  units <- vapply(map$target,function(t)if(startsWith(t,"d"))spec$scale[as.integer(substring(t,2L))]^.np_bias_derivative_order(spec,as.integer(substring(t,2L))) else 1,numeric(1))
  C <- sweep(C,1L,units,"/"); L <- sweep(L,1L,units,"/")
  factor <- sweep(C,2L,sqrt(spec$variance),"*")
  if(any(!is.finite(factor)))stop("nonfinite corrected influence rows",call.=FALSE)
  list(W=C,original=L,center=drop(C%*%spec$y),original_center=drop(L%*%spec$y),
    se=sqrt(rowSums(factor^2)),factor=factor,main.regularization=regularization)
}
