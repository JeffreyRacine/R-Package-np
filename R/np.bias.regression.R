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
  out <- .np_semihat_make_regbw_state(source, x, numeric(nrow(x)), bandwidth,
    kernel.choices = c("gaussian", "epanechnikov", "uniform", "beta"))
  if (!identical(out[["ckertype"]], source[["ckertype"]]) ||
      !isTRUE(out[["ckerorder"]] == source[["ckerorder"]]) ||
      !identical(out[["ckerbound"]], source[["ckerbound"]]))
    stop("regression bias preparation changed the requested kernel identity", call. = FALSE)
  out[["call"]] <- NULL
  out
}

.np_bias_regression_nn_zero <- function(bws, kernel = "epanechnikov") {
  isTRUE(bws[["type"]] %in% c("generalized_nn", "adaptive_nn")) &&
    isTRUE(bws[["ckertype"]] %in% kernel) &&
    isTRUE(bws[["ckerorder"]] == 2L) &&
    isTRUE(bws[["ckerbound"]] %in% c("none", "fixed")) &&
    (identical(bws[["regtype"]], "lc") ||
     (identical(bws[["regtype"]], "lp") &&
      identical(bws[["basis"]], "glp") &&
      is.numeric(bws[["degree"]]) && length(bws[["degree"]]) > 0L &&
      all(!is.na(bws[["degree"]]) & bws[["degree"]] == 0L)))
}

.np_bias_regression_nn_first_order <- function(bws, kernel) {
  isTRUE(bws[["type"]] %in% c("generalized_nn", "adaptive_nn")) &&
    isTRUE(bws[["ckertype"]] %in% kernel) &&
    isTRUE(bws[["ckerorder"]] == 2L) &&
    isTRUE(bws[["ckerbound"]] %in% c("none", "fixed")) &&
    (identical(bws[["regtype"]], "ll") ||
     (identical(bws[["regtype"]], "lp") &&
      identical(bws[["basis"]], "glp") &&
      is.numeric(bws[["degree"]]) && length(bws[["degree"]]) > 0L &&
      all(!is.na(bws[["degree"]]) & bws[["degree"]] %in% c(0L, 1L)) &&
      any(bws[["degree"]] == 1L)))
}

.np_bias_regression_kernel_class <- function(bws, calibration) {
  if (identical(bws[["ckertype"]], "gaussian") &&
      bws[["ckerorder"]] == 2L &&
      (identical(bws[["ckerbound"]], "none") ||
       (identical(bws[["type"]], "fixed") &&
        identical(bws[["ckerbound"]], "fixed") &&
        isTRUE(bws[["nuno"]] == 0L) && isTRUE(bws[["nord"]] == 0L)) ||
       .np_bias_regression_nn_zero(bws, "gaussian") ||
       .np_bias_regression_nn_first_order(bws, "gaussian")))
    return("incumbent")
  # Fixed constructions share the qualified complete non-Gaussian rows.
  # Degree-zero NN and first-order Epanechnikov NN share the qualified reducers.
  if ((identical(calibration, "basic-local-envelope") ||
       (identical(bws[["type"]], "fixed") &&
        isTRUE(calibration %in% c("basic-envelope", "separate-tails",
                                 "tail-envelope", "absolute"))) ||
       ((.np_bias_regression_nn_zero(bws, c("epanechnikov", "uniform", "beta")) ||
        .np_bias_regression_nn_first_order(bws, "epanechnikov")) &&
        isTRUE(calibration %in% c("basic-envelope", "separate-tails",
                                 "tail-envelope", "absolute")))) &&
      isTRUE(bws[["type"]] %in% c("fixed", "generalized_nn", "adaptive_nn")) &&
      bws[["ckerorder"]] == 2L &&
      ((identical(bws[["ckertype"]], "epanechnikov") &&
        (identical(bws[["ckerbound"]], "none") ||
         ((identical(bws[["type"]], "fixed") ||
           .np_bias_regression_nn_zero(bws) ||
           .np_bias_regression_nn_first_order(bws, "epanechnikov")) &&
          identical(bws[["ckerbound"]], "fixed")))) ||
       (identical(bws[["ckertype"]], "beta") &&
        (identical(bws[["type"]], "fixed") ||
         .np_bias_regression_nn_zero(bws, "beta")) &&
        identical(bws[["ckerbound"]], "fixed")) ||
       (identical(bws[["ckertype"]], "uniform") &&
        .np_bias_regression_nn_zero(bws, "uniform") &&
        isTRUE(bws[["ckerbound"]] %in% c("none", "fixed")))))
    return("local-envelope")
  "unsupported"
}

.np_bias_regression_calibration <- function(bws, calibration, method, band,
                                             bootstrap, wild, gradients = FALSE) {
  # These fixed fitted rows retain the scope of the previously supported
  # explicit local envelope until the physical-coordinate discrepancies are resolved.
  if (identical(bws[["type"]], "fixed") &&
      !identical(calibration, "basic-local-envelope")) {
    degree <- bws[["degree"]]
    raw.quadratic <- identical(bws[["ckertype"]], "epanechnikov") &&
      identical(bws[["regtype"]], "lp") &&
      identical(bws[["bernstein.basis"]], FALSE) && any(degree == 2L)
    beta.derivative <- identical(bws[["ckertype"]], "beta") && isTRUE(gradients) &&
      (identical(bws[["regtype"]], "ll") || any(degree > 0L))
    if (raw.quadratic || beta.derivative)
      stop("this fixed-bandwidth iterated regression target requires the explicit basic-local-envelope calibration for supported pairs or Mammen bands", call. = FALSE)
  }
  if (!identical(calibration, "default"))
    return(.np_plot_calibration_method(calibration, "regression"))
  # Resolve once from the fitted model and requested law, before preparing
  # the complete family. Explicit controls never delegate this decision.
  pairs <- bootstrap %in% c("inid", "fixed", "geom")
  mammen <- identical(bootstrap, "wild") && identical(wild, "mammen")
  if (identical(method, "bootstrap") && !identical(band, "pmzsd") &&
      (pairs || mammen) && isTRUE(bws[["ncon"]] > 0L) &&
      (.np_bias_regression_nn_zero(bws, c("epanechnikov", "uniform", "beta")) ||
        .np_bias_regression_nn_first_order(bws, "epanechnikov")))
    return("basic-local-envelope")
  if (identical(method, "bootstrap") && !identical(band, "pmzsd") &&
      identical(bws[["type", exact=TRUE]], "fixed") &&
      isTRUE(bws[["ncon", exact=TRUE]] > 0L) && (pairs || mammen) &&
      identical(.np_bias_regression_kernel_class(bws, "basic-local-envelope"),
                "local-envelope"))
    return("basic-local-envelope")
  # For the distinct Rademacher law, retain both its signed and absolute
  # intervals. Resolve this declared-law choice once for the complete family.
  if (identical(method, "bootstrap") && !identical(band, "pmzsd") &&
      (identical(bws[["type", exact=TRUE]], "fixed") ||
       (.np_bias_regression_nn_zero(bws, c("epanechnikov", "uniform", "beta")) ||
        .np_bias_regression_nn_first_order(bws, "epanechnikov"))) &&
      isTRUE(bws[["ncon", exact=TRUE]] > 0L) &&
      identical(bootstrap, "wild") && identical(wild, "rademacher") &&
      identical(.np_bias_regression_kernel_class(bws, "tail-envelope"),
                "local-envelope"))
    return("tail-envelope")
  "absolute"
}

.np_bias_regression_kernel_law <- function(admission, law, bws = NULL) {
  fixed <- !is.null(admission) && identical(admission[["type"]], "fixed")
  nn.zero <- !is.null(admission) && !fixed &&
    .np_bias_regression_nn_zero(bws, c("epanechnikov", "uniform", "beta"))
  nn.first <- !is.null(admission) && !fixed &&
    .np_bias_regression_nn_first_order(bws, "epanechnikov")
  if (!is.null(admission) && !fixed) {
    allowed <- if (nn.zero || nn.first) c("gaussian", "mammen", "rademacher", "inid", "fixed", "geom") else
      c("inid", "fixed", "geom")
    if (!is.character(law) || length(law) != 1L || is.na(law) || !law %in% allowed)
      stop(if (nn.zero)
        "degree-zero kernel NN inference requires Gaussian, Mammen, Rademacher, inid, fixed or geom resampling" else if (nn.first)
        "first-order Epanechnikov NN inference requires Gaussian, Mammen, Rademacher, inid, fixed or geom resampling" else
        "second-order kernel NN envelopes require inid, fixed or geom resampling", call.=FALSE)
  }
  if (fixed || nn.zero || nn.first) {
    description <- if (fixed) "fixed non-Gaussian" else if (nn.zero) "degree-zero NN non-Gaussian" else
      "first-order Epanechnikov NN"
    if (identical(admission[["calibration"]], "absolute") &&
        !isTRUE(law %in% c("gaussian", "mammen", "rademacher", "inid", "fixed", "geom")))
      stop(paste(description, "absolute calibration requires Gaussian, Mammen, Rademacher, inid, fixed or geom resampling"), call.=FALSE)
    if (isTRUE(admission[["calibration"]] %in% c("separate-tails", "tail-envelope")) &&
        !isTRUE(law %in% c("mammen", "rademacher", "inid", "fixed", "geom")))
      stop(paste(description, "tail calibration requires Mammen, Rademacher, inid, fixed or geom resampling"), call.=FALSE)
  }
  invisible(NULL)
}

# Bounds and data must use the same arithmetic at exact support endpoints.
.np_bias_regression_support <- function(bws, continuous, offset, scale) {
  if (identical(bws[["ckerbound"]], "fixed")) {
    bws[["ckerlb"]][continuous] <- (bws[["ckerlb"]][continuous] - offset) / scale
    bws[["ckerub"]][continuous] <- (bws[["ckerub"]][continuous] - offset) / scale
  }
  bws
}

.np_bias_regression_prepare <- function(bws, x, y, control = np_bias_control(), adapt = FALSE, gradient.order = NULL, law = NULL, calibration = NULL) {
  x <- toFrame(x)
  if (!is.numeric(y) || length(y) != nrow(x) || any(!is.finite(y)) || anyNA(x))
    stop("bias correction requires finite retained training data", call. = FALSE)
  reg <- npValidatedConditionalRegSpec(bws, where = "bias correction", ncon.field = "ncon")
  nc <- which(bws[["icon"]]); fc <- which(!bws[["icon"]]); d <- length(nc)
  if (!d) return(.np_bias_categorical_prepare(bws,x,y,control,reg))
  legacy.hint <- if (bws[["type"]]=="fixed" && d>0L)
    "; explicit bootstrap bias.control=np_bias_control(construction='legacy') retains the previous construction for supported bootstrap methods" else ""
  kernel.class <- .np_bias_regression_kernel_class(bws, calibration)
  if (kernel.class == "unsupported" || reg[["basis.engine"]] != "glp")
    stop(paste0("iterated regression bias correction requires a supported second-order kernel, at least one continuous predictor and a GLP main model; fixed-support Gaussian correction requires fixed bandwidths with continuous predictors only, or qualified nearest-neighbour models with main degrees at most one",legacy.hint), call. = FALSE)
  if (length(fc) && !all(vapply(x[fc],is.factor,logical(1))))
    stop("bias correction requires factor-valued categorical predictors",call.=FALSE)
  admission <- if (kernel.class == "local-envelope") list(
    version="regression-kernel-envelope-v1", calibration=calibration,
    kernel=bws[["ckertype"]], order=bws[["ckerorder"]], type=bws[["type"]],
    bound=bws[["ckerbound"]], lower=bws[["ckerlb"]], upper=bws[["ckerub"]]) else NULL
  .np_bias_regression_kernel_law(admission, law, bws)
  if (!is.null(admission) && identical(bws[["type"]], "generalized_nn") &&
      any(reg[["degree.engine"]] > 1L))
    stop("second-order kernel generalized-NN envelopes are currently qualified only for main degrees at most one", call.=FALSE)
  orders <- .np_bias_derivative_request(reg,d,gradient.order)
  if (.np_bias_profile_layout(bws,x)) {
    spec <- .np_bias_profile_prepare(bws,x,y,control,reg,orders,adapt,law)
    if (!is.null(admission)) spec$settings$kernel.admission <- admission
    return(spec)
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
  main.source <- .np_bias_regression_support(bws,nc,offset,scale)
  main <- .np_bias_regbw(main.source,standardized,h)
  spec <- list(x=standardized,xcon=standardized[nc],numeric.cols=nc,factor.cols=fc,
    d=d,q=q,bernstein=settings[["bernstein.basis"]],keys=keys,groups=groups,bws=main,
    b=unname(settings[["bandwidth"]][["p0"]]/scale),
    h=unname(settings[["bandwidth"]][["operator"]]/scale),
    variance_h=unname(settings[["bandwidth"]][["variance"]]/scale),
    offset=offset,scale=scale,y=as.double(y),settings=settings)
  if (!is.null(admission)) spec$settings$kernel.admission <- admission
  # Degree-zero uniform support ties use the fitted physical coordinates.
  # Auxiliary Gaussian projections retain their standardized coordinates.
  if (.np_bias_regression_nn_zero(bws, "uniform"))
    spec[["main.physical"]] <- list(x=x,
      bws=.np_bias_regbw(bws,x,unname(bws[["bandwidth"]][["x"]])))
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
    if (.np_bias_local_supported(spec,law)) spec <- .np_bias_local_prepare(spec) else
      spec$variance <- .np_bias_variance(spec,spec$residual)
  }
  if (!is.null(policy)) spec$settings$adaptation <- policy
  if (.np_bias_main_capture_supported(spec,law)) {
    if (law %in% c("inid","fixed","geom")) {
      spec$main.capture <- list(version="fixed-pairs-fitted-v1",law=law)
      spec$settings$main.capture <- list(version="fixed-pairs-fitted-v1",law=law,
        main="complete original fitted capture",variance="unchanged smoothed P0 residual variance",
        qualification="continuous fixed Gaussian2, inactive adaptation, mean/first derivative, inid/fixed/geom",
        storage="no fixed capture ceiling; bounded tiles")
    } else {
    spec$main.capture <- list(version="fixed-multivariate-fitted-v1")
    spec$settings$main.capture <- list(version="fixed-multivariate-fitted-v1",
      main="complete original fitted capture",variance="unchanged smoothed P0 residual variance",
      qualification="continuous fixed Gaussian2, inactive adaptation, mean/first derivative, Gaussian/Mammen")
    }
  }
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
  .np_bias_main_admit(spec,H,ridge,reg)
}

.np_bias_main_admit <- function(spec,H,ridge,reg,captured=FALSE) {
  if (is.null(ridge) || any(!is.finite(ridge)) || any(ridge<0) || any(!is.finite(H)))
    stop("bias correction requires finite canonical main observation rows",call.=FALSE)
  admitted <- spec$bws[["type"]]=="generalized_nn" &&
    any(reg[["degree.engine"]]>1L) && isTRUE(reg[["bernstein.basis.engine"]])
  admitted <- admitted || (isTRUE(captured) && !is.null(spec$local.inference))
  admitted <- admitted || (!is.null(spec$main.capture) && inherits(captured,"np_bias_main_capture_certificate") &&
    identical(captured$state,spec))
  if (any(ridge!=0) && !admitted)
    stop("regularized main rows require a qualified fitted capture or generalized-NN GLP model with a stored Bernstein basis; the fitted model has not been changed",call.=FALSE)
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

.np_bias_regression_rows <- function(spec,grid,map,original=NULL,captured=FALSE) {
  # This admission covers means and differences of means. Keep continuous
  # uniform derivatives outside it, including the pooled dispatch below.
  if (identical(spec[["settings"]][["main"]][["kernel"]], "uniform") &&
      any(startsWith(map[["target"]], "d")))
    stop("iterated uniform-kernel correction is qualified only for means and categorical mean contrasts; continuous derivatives remain unqualified; explicit bias.control=np_bias_control(construction='legacy') retains previous supported bootstrap routes", call.=FALSE)
  if (identical(spec[["settings"]][["main"]][["kernel"]], "beta") &&
      isTRUE(spec[["settings"]][["main"]][["type"]] %in%
             c("generalized_nn", "adaptive_nn")) &&
      any(startsWith(map[["target"]], "d")))
    stop("iterated degree-zero beta-kernel NN correction is qualified only for means and categorical mean contrasts; continuous derivatives remain unqualified", call.=FALSE)
  if (!is.null(original) && (isTRUE(spec$categorical.only) || isTRUE(spec$pooled)))
    stop("captured original rows are unavailable for this bias owner",call.=FALSE)
  if (isTRUE(spec$categorical.only)) return(.np_bias_categorical_rows(spec,grid,map))
  if (isTRUE(spec$pooled)) return(.np_bias_profile_rows(spec,grid,map))
  if (nrow(map)>500L || nrow(unique(grid[spec$numeric.cols]))>25L)
    stop("internal bias target block exceeded",call.=FALSE)
  physical <- spec[["main.physical"]]
  main.grid <- if (!is.null(physical)) grid else NULL
  for(j in seq_len(spec$d))grid[[spec$numeric.cols[j]]] <- (grid[[spec$numeric.cols[j]]]-spec$offset[j])/spec$scale[j]
  L <- if (is.null(original)) {
    if (is.null(physical)) .np_bias_original(spec,grid,map) else {
      main <- spec
      main$x <- physical[["x"]]
      main$xcon <- main$x[spec$numeric.cols]
      main$bws <- physical[["bws"]]
      .np_bias_original(main,main.grid,map)
    }
  } else {
    diagnostic <- attr(original,"main.regularization",exact=TRUE)
    if (!is.matrix(original) || !is.numeric(original) ||
        !identical(dim(original),c(nrow(map),nrow(spec$x))) ||
        !is.data.frame(diagnostic) || nrow(diagnostic)!=nrow(map))
      stop("invalid bounded captured original rows",call.=FALSE)
    reg <- npValidatedConditionalRegSpec(spec$bws,where="bias main operator",ncon.field="ncon")
    .np_bias_main_admit(spec,original,diagnostic$ridge,reg,captured=captured)
    original
  }
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
