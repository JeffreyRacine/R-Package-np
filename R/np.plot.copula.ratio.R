# New simultaneous requests only. Compose existing NN numerator/marginal
# factors at the unchanged inverse coordinates; no kernel or radius owner here.
.np_plot_copula_ratio_positive <- function(a, b, operation = c("add", "multiply", "divide")) {
  operation <- match.arg(operation)
  if (any(!is.finite(a)) || any(!is.finite(b)) || any(a < 0) || any(b < 0) ||
      (identical(operation, "divide") && any(b <= 0)))
    stop("invalid copula ratio discrepancy operands", call. = FALSE)
  value <- switch(operation, add = a + b, multiply = a * b, divide = a / b)
  lost <- switch(operation, add = (a > 0 | b > 0) & value == 0,
    multiply = a > 0 & b > 0 & value == 0, divide = a > 0 & value == 0)
  if (any(!is.finite(value)) || any(lost))
    .np_plot_joint_unavailable("copula ratio discrepancy is not representable")
  value
}

.np_plot_copula_ratio_combine <- function(numerator, marginal, point,
                                           numerator.error, marginal.error,
                                           point.error, identity = FALSE) {
  if (!is.matrix(numerator) || !identical(dim(numerator), dim(marginal)) ||
      !identical(dim(numerator), dim(numerator.error)) ||
      !identical(dim(numerator), dim(marginal.error)) ||
      length(point) != nrow(numerator) || length(point.error) != nrow(numerator) ||
      any(!is.finite(c(numerator, marginal, point, numerator.error,
                      marginal.error, point.error))) ||
      any(numerator.error < 0) || any(marginal.error < 0) || any(point.error < 0))
    stop("invalid copula ratio composition", call. = FALSE)
  if (isTRUE(identity)) {
    # Caller has reused the SAME one-dimensional numerator primitive as its
    # sole marginal. Numerical equality alone is not this certificate.
    if (!identical(numerator, marginal) || any(point != 1))
      stop("one-dimensional copula identity contract failed", call. = FALSE)
    return(list(A = numerator * 0, error = numerator * 0))
  }
  add <- function(a,b) .np_plot_copula_ratio_positive(a,b,"add")
  mul <- function(a,b) .np_plot_copula_ratio_positive(a,b,"multiply")
  correction <- point * marginal
  A <- numerator - correction
  if (any(!is.finite(c(correction,A))))
    .np_plot_joint_unavailable("copula ratio factor is not representable")
  error <- add(numerator.error, add(mul(abs(point),marginal.error),
    mul(point.error,add(abs(marginal),marginal.error))))
  # Arithmetic allowance is relative to the actual summands, never to one.
  error <- add(error, mul(8*.Machine$double.eps,
    add(abs(numerator),abs(correction))))
  list(A=A,error=error)
}

.np_plot_copula_ratio_influence <- function(object, data, ex, se,
                                            owner = environment(.np_plot_copula_ratio_influence),
                                            max.bytes = 128*1024^2,
                                            tile.rows = NULL) {
  if (!is.environment(owner)) stop("copula ratio owner must be an environment")
  f <- function(name) get(name,envir=owner,mode="function",inherits=TRUE)
  unavailable <- f(".np_plot_joint_unavailable")
  bws <- object[["bws",exact=TRUE]]
  type <- bws[["type",exact=TRUE]]
  if (!identical(object[["target",exact=TRUE]],"density") ||
      !type %in% c("generalized_nn","adaptive_nn") ||
      !is.data.frame(data) || !is.data.frame(ex) || anyNA(data) || anyNA(ex) ||
      !identical(names(data),bws[["xnames",exact=TRUE]]) ||
      !identical(names(ex),names(data)))
    stop("invalid NN copula ratio data contract",call.=FALSE)
  n <- nrow(data); m <- nrow(ex); p <- ncol(data)
  incumbent <- object[["copula",exact=TRUE]]
  if (n<2L || m<1L || p<1L || !is.numeric(incumbent) || length(incumbent)!=m ||
      any(!is.finite(incumbent)) || length(max.bytes)!=1L ||
      !is.finite(max.bytes) || max.bytes<=0)
    stop("invalid NN copula ratio dimensions or budget",call.=FALSE)
  tile <- min(m,64L,max(1L,floor(4*1024^2/(8*as.double(n)*(80+12*p)))))
  if (!is.null(tile.rows)) {
    if (length(tile.rows)!=1L || !is.finite(tile.rows) || tile.rows!=floor(tile.rows) ||
        tile.rows<1L || tile.rows>tile)
      stop("tile.rows may only reduce the positive copula ratio tile",call.=FALSE)
    tile <- as.integer(tile.rows)
  }
  # Final factors/errors and calibration copies, a held numerator tile,
  # accumulated margins, the current child workspace and O(n*p) metadata.
  parent.bytes <- 48*as.double(n)*m + 8*as.double(n)*tile*(24+2*p) +
    512*as.double(n)*(p+1) + 256*as.double(m)*(p+1)
  child.bytes <- 64*1024^2 + 48*as.double(n)*tile +
    8*as.double(n)*tile*(20+8*p) + 256*as.double(n)*(p+1) +
    128*as.double(tile)*(p+1)
  forecast <- parent.bytes+child.bytes
  if (!is.finite(forecast) || forecast>max.bytes)
    unavailable("NN copula ratio exceeds the aggregate plot memory budget")
  child.limit <- max.bytes-parent.bytes
  ann <- identical(type,"adaptive_nn") && any(bws[["icon",exact=TRUE]])
  add <- function(a,b) f(".np_plot_copula_ratio_positive")(a,b,"add")
  mul <- function(a,b) f(".np_plot_copula_ratio_positive")(a,b,"multiply")
  div <- function(a,b) f(".np_plot_copula_ratio_positive")(a,b,"divide")
  normal <- f(".np_plot_density_unconditional")
  ann.owner <- f(".np_plot_ann_bounded_density_factor")
  mb <- lapply(seq_len(p),function(j)
    f(".npcopula_marginal_bw")(bws,data,j,target="density",kbandwidth=FALSE))
  train <- lapply(seq_len(p),function(j)f(".npcopula_marginal_data")(bws,data,j))
  evaluate <- lapply(seq_len(p),function(j)f(".npcopula_marginal_eval_data")(bws,ex,j))
  if (p==1L &&
      (!identical(data,train[[1L]]) || !identical(ex,evaluate[[1L]]) ||
       !identical(f(".np_make_kbandwidth_unconditional")(bws,data),
                  f(".np_make_kbandwidth_unconditional")(mb[[1L]],train[[1L]]))))
    stop("one-dimensional copula marginal metadata changed the canonical primitive",call.=FALSE)
  A <- EA <- matrix(NA_real_,m,n)
  point <- point.error <- denominator <- numeric(m)
  progress <- f(".np_plot_bootstrap_progress_begin")(m,
    "Computing copula joint numerator and marginal influence",unit="row")
  on.exit(f(".np_plot_progress_end")(progress),add=TRUE)
  for (first in seq.int(1L,m,by=tile)) {
    ii <- seq.int(first,min(m,first+tile-1L)); nr <- length(ii)
    ee <- ex[ii,,drop=FALSE]
    # Convert only known numerical/regularity conditions. Malformed native
    # output, wrong owner or implementation failures remain hard errors.
    invoke.ann <- function(bb,dd,xx) tryCatch(
      ann.owner(bb,dd,xx,owner=owner,max.bytes=child.limit,tile.rows=nrow(xx)),
      error=function(e) {
        prefixes <- c("np_plot_ann_unavailable: ","np_plot_beta_numerical_unavailable: ",
                      "np_plot_cdf_numerical_unavailable: ","np_plot_nn_unavailable: ")
        found <- which(vapply(prefixes,function(x)startsWith(conditionMessage(e),x),logical(1L)))
        if(length(found)) unavailable(substring(conditionMessage(e),nchar(prefixes[found[1L]])+1L))
        stop(e)
      })
    joint <- if(ann) invoke.ann(bws,data,ee) else normal(bws,data,ee,cdf=FALSE)
    if (!is.matrix(joint$A) || !identical(dim(joint$A),c(nr,n)) ||
        any(!is.finite(joint$A)) || length(joint$estimate)!=nr ||
        any(!is.finite(joint$estimate)))
      stop("malformed NN copula numerator factor",call.=FALSE)
    # A one-dimensional copula density is a ratio of the same canonical
    # marginal primitive. Reuse it constructively, rather than compare two
    # approximate factors and guess that their difference is zero.
    if (p==1L) {
      if (any(incumbent[ii]!=1))
        unavailable("one-dimensional copula density does not certify its constant ratio")
      if (ann) {
        if (any(!is.finite(joint$denominator)) || any(joint$denominator<=0))
          unavailable("copula ratio requires positive raw marginal densities")
        denominator[ii] <- joint$denominator
        B0 <- joint$A
      } else {
        if(any(joint$estimate<=0))unavailable("copula ratio requires positive raw marginal densities")
        denominator[ii] <- joint$estimate
        B0 <- joint$A/joint$estimate
      }
      combined <- f(".np_plot_copula_ratio_combine")(B0,B0,rep(1,nr),
        matrix(0,nr,n),matrix(0,nr,n),numeric(nr),identity=TRUE)
      A[ii,] <- combined$A;EA[ii,] <- combined$error;point[ii] <- 1
      progress <- f(".np_plot_progress_tick")(progress,max(ii))
      next
    }
    total <- total.error <- matrix(0,nr,n)
    D <- rep(1,nr)
    for (j in seq_len(p)) {
      ej <- evaluate[[j]][ii,,drop=FALSE]
      if (ann && any(mb[[j]][["icon",exact=TRUE]])) {
        marginal <- invoke.ann(mb[[j]],train[[j]],ej)
        Bj <- marginal$A;Ej <- marginal$factor.absolute.discrepancy
        mj <- marginal$denominator
      } else {
        marginal <- normal(mb[[j]],train[[j]],ej,cdf=FALSE)
        mj <- marginal$estimate
        if(any(!is.finite(mj))||any(mj<=0))
          unavailable("copula ratio requires positive raw marginal densities")
        Bj <- marginal$A/mj
        if(ann)Bj <- Bj*sqrt(n/(n-1))
        Ej <- matrix(0,nr,n)
      }
      if (!is.matrix(Bj) || !identical(dim(Bj),c(nr,n)) ||
          !is.matrix(Ej) || !identical(dim(Ej),c(nr,n)) ||
          length(mj)!=nr || any(!is.finite(c(Bj,Ej,mj))) ||
          any(Ej<0) || any(mj<=0))
        stop("malformed NN copula marginal factor",call.=FALSE)
      D <- D*mj
      if(any(!is.finite(D))||any(D<=0))
        unavailable("NN copula marginal product is not representable")
      before <- total
      total <- total+Bj
      total.error <- add(total.error,add(Ej,
        mul(4*.Machine$double.eps,add(abs(before),abs(Bj)))))
      marginal <- Bj <- Ej <- before <- NULL
    }
    if(ann) {
      if(length(joint$denominator)!=nr ||
          any(abs(D-joint$denominator)>2e-11*abs(joint$denominator)))
        stop("ANN marginal factors changed the canonical denominator",call.=FALSE)
      B0 <- joint$A;E0 <- joint$factor.absolute.discrepancy
      cvalue <- joint$estimate;ec <- joint$point.absolute.discrepancy
      if(is.null(ec)||length(ec)!=nr||any(!is.finite(ec))||any(ec<0))
        stop("ANN ratio requires its point-discrepancy contract",call.=FALSE)
    } else {
      B0 <- joint$A/D;E0 <- matrix(0,nr,n)
      cvalue <- joint$estimate/D;ec <- numeric(nr)
    }
    combined <- f(".np_plot_copula_ratio_combine")(B0,total,cvalue,E0,total.error,ec)
    A[ii,] <- combined$A;EA[ii,] <- combined$error
    point[ii] <- cvalue;point.error[ii] <- ec;denominator[ii] <- D
    joint <- combined <- B0 <- E0 <- total <- total.error <- NULL
    progress <- f(".np_plot_progress_tick")(progress,max(ii))
  }
  scale <- f(".np_plot_four_factor_se")(A)
  error <- f(".np_plot_four_factor_se")(EA)
  relative.error <- ifelse(scale>0,error/scale,ifelse(error==0,0,Inf))
  if(any(!is.finite(relative.error))||any(relative.error>1e-8))
    unavailable("NN copula ratio numerical discrepancy exceeds its factor-relative budget")
  roundoff <- max(128,8*as.double(n))*.Machine$double.eps*pmax(abs(point),abs(incumbent))
  if(any(!is.finite(point)) || any(abs(point-incumbent)>.01*scale+roundoff))
    stop("NN copula ratio does not reproduce the plotted point",call.=FALSE)
  list(A=A,estimate=point,simultaneous.se=scale,exact.diagonal=TRUE,
    denominator=denominator,point.absolute.discrepancy=point.error,
    factor.absolute.discrepancy=EA,factor.discrepancy.norm=error,
    factor.discrepancy.ratio=relative.error,forecast.bytes=forecast,tile.rows=tile,
    order=list(evaluation=seq_len(m),observation=seq_len(n)),
    target=paste("NN copula density; full same-sample numerator/marginal ratio at frozen inverse coordinates;",
      if(ann)"ANN random-radius working influence" else "fixed or frozen-query-radius empirical-mass influence"),
    convention=list(inverse.coordinates="frozen",denominator="same-sample estimated marginals",
      geometry=if(ann)"same-sample ANN random-radius approximation" else "fixed categorical weights or frozen GNN query radii"),
    marginal.contract="new simultaneous factor scale; incumbent pointwise and Bonferroni intervals retained",
    numerical.contract="inherited primitive accuracy plus same-shape ANN PDF agreement and final ratio propagation; no universal special-function enclosure")
}
