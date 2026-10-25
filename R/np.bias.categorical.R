# Exact finite-profile reduction of the response-linear correction when q=0.
# The saturated profile projection annihilates every residual update; neither
# continuous pilot geometry nor a training-square projection is materialized.
.np_bias_categorical_prepare <- function(bws, x, y, control, reg) {
  control <- .np_bias_validate_control(control)
  if (control[["construction"]] != "iterated" || bws[["type"]] != "fixed" ||
      !ncol(x) || !all(vapply(x, is.factor, logical(1))))
    stop("categorical-only bias correction requires fixed-bandwidth factor predictors", call.=FALSE)
  geometry <- c("degree", "bernstein.basis",
    paste0(rep(c("p0", "p1", "operator", "variance"), each=2L), c(".scale", ".bandwidth")))
  if (length(intersect(control[["supplied"]], geometry)) || !is.null(control[["degree"]]) ||
      !identical(control[["bandwidth.controls"]], np_bias_control()[["bandwidth.controls"]]))
    stop("continuous auxiliary geometry controls do not apply to categorical-only bias correction", call.=FALSE)
  if (identical(control[["adaptation"]], "local-to-global"))
    stop("local-to-global adaptation requires a continuous predictor", call.=FALSE)
  fc <- seq_len(ncol(x)); keys <- .np_bias_profile_keys(x, fc)
  groups <- split(seq_len(nrow(x)), keys); sizes <- lengths(groups)
  if (any(sizes<2L))
    stop("categorical-only bias correction requires residual support in every observed profile", call.=FALSE)
  means <- vapply(groups, function(ii) mean(y[ii]), numeric(1))
  variances <- vapply(groups, function(ii) var(y[ii]), numeric(1))
  if (any(!is.finite(variances)) || any(variances<=0))
    stop("categorical-only bias correction requires positive within-profile variance", call.=FALSE)
  index <- match(keys, names(groups)); mu <- unname(means[index])
  q <- unname(1-1/sizes[index]); raw <- y-mu
  settings <- list(construction="iterated", version=3L, policy="finite-profile-v1",
    requested=control, iterations=control[["iterations"]], columns=character(),
    n=nrow(x), profile_sizes=sizes, n_min=min(sizes),
    auxiliary_degree=integer(), bandwidth=list(), effective_package_scale=list(),
    bernstein.basis=control[["bernstein.basis"]],
    auxiliary_space=list(construction="observed profile indicators", profiles=names(groups),
      columns=length(groups), materialized=FALSE),
    iteration_reduction="all residual increments are algebraically zero",
    studentization=list(gaussian="profile variance Gaussian errors",
      mammen="leverage-normalized residual wild errors"),
    adaptation=list(requested=control[["adaptation"]],active=FALSE,eligible=FALSE,
      reason="no continuous predictors"),
    training=list(n=nrow(x),columns=names(x),row.names=rownames(x),
      rows.omit=bws[["rows.omit"]],source="retained plot training frame"),
    main=list(bandwidth=bws[["bandwidth"]][["x"]],type=bws[["type"]],
      regtype=bws[["regtype"]],degree=bws[["degree"]],basis=bws[["basis"]],
      bernstein.basis=bws[["bernstein.basis"]],kernel=bws[["ckertype"]],
      order=bws[["ckerorder"]],scale.factor=bws[["sfactor"]][["x"]]))
  list(categorical.only=TRUE,x=x,xcon=x[,FALSE,drop=FALSE],y=as.double(y),bws=bws,
    d=0L,numeric.cols=integer(),factor.cols=fc,keys=keys,groups=groups,
    profile.index=index,mu=mu,residual_raw=raw,residual_row_variance=q,
    residual=raw/sqrt(q),variance=unname(variances[index]),variance_h=numeric(),settings=settings)
}

.np_bias_categorical_average_rows <- function(groups, index, n) {
  if (length(index)>500L || anyNA(index) || any(index<1L | index>length(groups)))
    stop("invalid bounded categorical profile row request", call.=FALSE)
  out <- matrix(0,length(index),n)
  for (j in seq_along(index)) {
    ii <- groups[[index[j]]]; out[j,ii] <- 1/length(ii)
  }
  out
}

.np_bias_categorical_student_prepare <- function(state) {
  list(pilot=list(x=state$xcon,groups=state$groups,categorical.only=TRUE,
      profile.index=state$profile.index),q=state$residual_row_variance,
    rmu=rep(0,nrow(state$x)),root.variance=sqrt(state$variance),
    wild.residual=state$residual,
    smoother=list(x=state$xcon,groups=state$groups,h=numeric()))
}

.np_bias_categorical_rows <- function(state, grid, map) {
  if (nrow(map)>500L) stop("categorical functional block exceeded", call.=FALSE)
  frames <- lapply(seq_len(nrow(map)),function(r)
    .np_bias_frames(state,grid[map$row[r],,drop=FALSE],map$target[r]))
  at <- do.call(rbind,lapply(frames,function(f)do.call(rbind,f$e)))
  index <- match(.np_bias_profile_keys(at,state$factor.cols),names(state$groups))
  if (anyNA(index)) stop("categorical-only bias correction requires observed endpoint profiles",call.=FALSE)
  owner <- rep(seq_along(frames),lengths(lapply(frames,`[[`,"e")))
  signs <- unlist(lapply(frames,`[[`,"sign"))
  C <- L <- matrix(0,nrow(map),nrow(state$x)); diagnostics <- list()
  for (ii in .np_bias_profile_tiles(nrow(at),25L)) {
    H <- .np_bias_categorical_average_rows(state$groups,index[ii],nrow(state$x))
    M <- npreghat(state$bws,txdat=state$x,exdat=at[ii,,drop=FALSE],
      output="matrix",.np.require.finite=TRUE)
    ridge <- attr(M,"ridge.used",exact=TRUE)
    if (length(ridge)!=length(ii) || any(!is.finite(ridge)) || any(ridge!=0))
      stop("categorical-only correction requires unregularized main rows",call.=FALSE)
    for (j in seq_along(ii)) {
      k <- owner[ii[j]]; C[k,] <- C[k,]+signs[ii[j]]*H[j,]
      L[k,] <- L[k,]+signs[ii[j]]*M[j,]
    }
    diagnostics[[length(diagnostics)+1L]] <- data.frame(row=map$row[owner[ii]],
      target=map$target[owner[ii]],frame=sequence(lengths(lapply(frames,`[[`,"e")))[ii],ridge=ridge)
  }
  F <- sweep(C,2L,sqrt(state$variance),"*"); se <- sqrt(rowSums(F^2))
  if (any(!is.finite(F))) stop("nonfinite categorical corrected rows",call.=FALSE)
  list(W=C,original=L,center=drop(C%*%state$y),original_center=drop(L%*%state$y),
    se=se,factor=F,main.regularization=do.call(rbind,diagnostics))
}
